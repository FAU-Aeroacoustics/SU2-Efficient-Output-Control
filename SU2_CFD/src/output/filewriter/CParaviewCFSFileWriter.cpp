/*!
 * \file CParaviewCFSFileWriter.cpp
 * \brief Filewriter class for Paraview binary format.
 * \author E. Bagheri
 * \version 7.3.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../../../include/output/filewriter/CParaviewCFSFileWriter.hpp"
#include "../../../../Common/include/toolboxes/printing_toolbox.hpp"
#include "../../../include/output/highfive/H5Easy.hpp"
#include "../../../include/output/highfive/H5File.hpp"

const string CParaviewCFSFileWriter::fileExt = ".h5";
bool CParaviewCFSFileWriter::topoState = false;
int CParaviewCFSFileWriter::steadyIteration = 0;

CParaviewCFSFileWriter::CParaviewCFSFileWriter(CParallelDataSorter* valDataSorter, CConfig* config)
    : TimeIter(config->GetTimeIter()),
      lastTimeStep(config->GetnTime_Iter()),
      timeStep(SU2_TYPE::GetValue(config->GetTime_Step())),
      CFileWriter(valDataSorter, fileExt) {
  if (config->GetTime_Marching() == TIME_MARCHING::STEADY) {
    steadySim = true;
  }

  else {
    steadySim = false;
  }

  if (rank == MASTER_NODE && isSteady()) {
    incrementSteadyIteration();
  }
}

CParaviewCFSFileWriter::~CParaviewCFSFileWriter() {}

void CParaviewCFSFileWriter::incrementSteadyIteration() { steadyIteration++; };

void CParaviewCFSFileWriter::Write_Data(string val_filename) {
  if (!dataSorter->GetConnectivitySorted()) {
    SU2_MPI::Error("Connectivity must be sorted.", CURRENT_FUNCTION);
  }
  unsigned short lastindex = val_filename.find_last_of("_");

  std::string H5_fileName;
  std::string timestep;

  /*--- create the name for the output file ---*/
  if (isSteady()) {
    timestep = "Step_0";
    H5_fileName = val_filename + CParaviewCFSFileWriter::fileExt;
  } else {
    timestep = "Step_" + std::to_string(TimeIter);
    H5_fileName = val_filename.substr(0, lastindex) + CParaviewCFSFileWriter::fileExt;
  }

  const int NCOORDS = 3;
  const unsigned short nDim = dataSorter->GetnDim();
  unsigned short iDim = 0;

  /*--- Array containing the field names we want to output ---*/

  const vector<string>& fieldNames = dataSorter->GetFieldNames();

  unsigned long iPoint, iElem;
  unsigned long myPoint, GlobalPoint;

  GlobalPoint = dataSorter->GetnPointsGlobal();
  myPoint = dataSorter->GetnPoints();

  /*--- Compute our local number of elements, the required storage,
   and reduce the total number of elements and storage globally. ---*/

  unsigned long myElem, myElemStorage, GlobalElem, GlobalElemStorage;

  unsigned long nParallel_Line = dataSorter->GetnElem(LINE), nParallel_Tria = dataSorter->GetnElem(TRIANGLE),
                nParallel_Quad = dataSorter->GetnElem(QUADRILATERAL),
                nParallel_Tetr = dataSorter->GetnElem(TETRAHEDRON), nParallel_Hexa = dataSorter->GetnElem(HEXAHEDRON),
                nParallel_Pris = dataSorter->GetnElem(PRISM), nParallel_Pyra = dataSorter->GetnElem(PYRAMID);

  myElem = dataSorter->GetnElem();
  myElemStorage = dataSorter->GetnConn();
  GlobalElem = dataSorter->GetnElemGlobal();
  GlobalElemStorage = dataSorter->GetnConnGlobal();

  unsigned short varStart = 2;
  if (nDim == 3) varStart++;

  /*--- Loop over all variables that have been registered in the output. ---*/
  unsigned short iField, VarCounter = varStart;
  for (iField = varStart; iField < fieldNames.size(); iField++) {
    string fieldname = fieldNames[iField];
    fieldname.erase(remove(fieldname.begin(), fieldname.end(), '"'), fieldname.end());

    /*--- Check whether this field is a vector or scalar. ---*/
    bool output_variable = true, isVector = false;
    size_t found = fieldNames[iField].find("_x");
    if (found != string::npos) {
      output_variable = true;
      isVector = true;
    }
    found = fieldNames[iField].find("_y");
    if (found != string::npos) {
      /*--- We have found a vector, so skip the Y component. ---*/
      output_variable = false;
      VarCounter++;
    }
    found = fieldNames[iField].find("_z");
    if (found != string::npos) {
      /*--- We have found a vector, so skip the Z component. ---*/
      output_variable = false;
      VarCounter++;
    }

    /*--- Write the point data as an <X,Y,Z> vector or a scalar. ---*/

    if (output_variable && isVector) {
      /*--- Adjust the string name to remove the leading "X-" ---*/

      fieldname.erase(fieldname.end() - 2, fieldname.end());

    } else if (output_variable) {
    }
  }

  vector<float> dataBufferFloat(myPoint * NCOORDS);
  vector<vector<double> > twoDimBuffer(myPoint, vector<double>(NCOORDS, 0.0));

  for (iPoint = 0; iPoint < myPoint; iPoint++) {
    for (iDim = 0; iDim < NCOORDS; iDim++) {
      if (nDim == 2 && iDim == 2) {
        twoDimBuffer[iPoint][iDim] = 0.0;
        dataBufferFloat[iPoint * NCOORDS + iDim] = 0.0;
      } else {
        float val = (float)dataSorter->GetData(iDim, iPoint);
        twoDimBuffer[iPoint][iDim] = val;
        dataBufferFloat[iPoint * NCOORDS + iDim] = val;
      }
    }
  }

  string MeshPathName = "/Mesh/Nodes/";
  string fullMeshPathName = MeshPathName + timestep + "/Coordinates";

  /*--- store the connectivity in a vector ---*/
  vector<vector<int> > connectivity(myElem, vector<int>(N_POINTS_HEXAHEDRON, 0));

  /*--- enum for CFS element types ---*/
  enum H5_Element_TYPE {
    H5VERTEX = 1,        /*!< \brief CFS element type for defining a vertex element. */
    H5LINE = 2,          /*!< \brief CFS element type for defining a line element. */
    H5TRIANGLE = 4,      /*!< \brief CFS element type for defining a triangle element. */
    H5QUADRILATERAL = 6, /*!< \brief CFS element type for defining a quadrilateral element. */
    H5TETRAHEDRON = 9,   /*!< \brief CFS element type for defining a tetrahedron element. */
    H5HEXAHEDRON = 11,   /*!< \brief CFS element type for defining a hexahedron element. */
    H5PRISM = 16,        /*!< \brief CFS element type for defining a prism element. */
    H5PYRAMID = 14       /*!< \brief CFS element type for defining a pyramid element. */
  };

  vector<int> connBuf(myElemStorage);
  vector<int> offsetBuf(myElem);
  vector<int> elementType(myElem, 0);

  unsigned long iStorage = 0, iElemID = 0;
  unsigned short iNode = 0;

  /*---write the mesh to the HDF5 ---*/
  if (!getTopoState()) {
    /*---attributes and result descriptions for paraview CFS file reader are written by the master node ---*/
    if (rank == MASTER_NODE) {
      HighFive::File file(H5_fileName,
                          HighFive::File::ReadWrite | HighFive::File::OpenOrCreate | HighFive::File::Truncate);
      string AnalysisType = "transient";

      HighFive::Group meshGroup = file.createGroup("Mesh");
      HighFive::Group meshElement = file.createGroup("Mesh/Elements");
      HighFive::Group meshNodes = file.createGroup("Mesh/Nodes");
      HighFive::Group meshRegions = file.createGroup("Mesh/Regions");
      HighFive::Group defaultRegions = file.createGroup("Mesh/Regions/default");
      HighFive::Group ResultsGrp = file.createGroup("Results");
      HighFive::Group ResultsMeshGrp = file.createGroup("Results/Mesh");
      HighFive::Group ResultsMeshMultstepGrp = file.createGroup("Results/Mesh/MultiStep_1");
      HighFive::Group ResultDescription = file.createGroup("Results/Mesh/MultiStep_1/ResultDescription");

      HighFive::Attribute resultMeshAttrib =
          ResultsMeshGrp.createAttribute<int>("ExternalFiles", HighFive::DataSpace(1));
      resultMeshAttrib.write(0);

      HighFive::Attribute multiStepAttrib1 =
          ResultsMeshMultstepGrp.createAttribute<std::string>("AnalysisType", HighFive::DataSpace::From(AnalysisType));
      multiStepAttrib1.write(AnalysisType);

      HighFive::Attribute multiStepAttrib2 =
          ResultsMeshMultstepGrp.createAttribute<int>("LastStepNum", HighFive::DataSpace(1));
      multiStepAttrib2.write(0);

      HighFive::Attribute defaultRegionsAttrib =
          defaultRegions.createAttribute<int>("Dimension", HighFive::DataSpace(1));
      defaultRegionsAttrib.write(3);
      std::vector<int> regionElems(GlobalElem);
      std::iota(std::begin(regionElems), std::end(regionElems), 1);  // Fill with 0, 1, ..., 99.

      std::vector<int> regionNodeList(GlobalPoint);
      std::iota(std::begin(regionNodeList), std::end(regionNodeList), 1);  // Fill with 0, 1, ..., 99.

      HighFive::Attribute meshElementAt1 = meshElement.createAttribute<int>("QuadraticElems", HighFive::DataSpace(1));
      meshElementAt1.write(0);

      meshElement.createAttribute<int>("Num_QUAD4", HighFive::DataSpace(1))
          .write(dataSorter->GetnElemGlobal(QUADRILATERAL));

      HighFive::Attribute meshAttrib = meshGroup.createAttribute<int>("Dimension", HighFive::DataSpace(1));
      meshAttrib.write(3);

      HighFive::DataSet dataSetRegionsEl =
          defaultRegions.createDataSet<int>("Elements", HighFive::DataSpace(GlobalElem));
      dataSetRegionsEl.select({0}, {GlobalElem}).write(regionElems);

      HighFive::DataSet dataSetRegionsNodes =
          defaultRegions.createDataSet<int>("Nodes", HighFive::DataSpace(GlobalPoint));
      dataSetRegionsNodes.select({0}, {GlobalPoint}).write(regionNodeList);

      unsigned short varStart = 2;
      if (nDim == 3) varStart++;

      for (iField = varStart; iField < fieldNames.size(); iField++) {
        std::vector<std::string> DOFNames = {""};
        std::vector<std::string> EntityName = {"default"};
        std::vector<std::string> unit = {"unknown"};
        std::vector<double> stepValues(getLastTimeStep(), 0);
        std::vector<int> stepNumbers(getLastTimeStep(), 0);
        std::iota(std::begin(stepNumbers), std::end(stepNumbers), 0);  // Fill with 0, 1, ..., 99.
        int entryType = 1;
        bool writeFieldName = true;
        passivedouble outputTimeStep = getTime_Step();
        string thisFieldName = fieldNames[iField];
        thisFieldName.erase(remove(thisFieldName.begin(), thisFieldName.end(), '"'), thisFieldName.end());

        for (unsigned int iStep = 0; iStep < stepValues.size(); iStep++) {
          stepValues[iStep] = stepNumbers[iStep] * outputTimeStep;
        }

        /*--- Check whether this field is a vector or scalar. ---*/
        bool output_variable = true, isVector = false;
        size_t found = fieldNames[iField].find("_x");
        if (found != string::npos) {
          entryType = 3;
          DOFNames[0] = "X";
          DOFNames.push_back("Y");
          DOFNames.push_back("Z");

          thisFieldName.erase(thisFieldName.end() - 2, thisFieldName.end());
        }

        found = fieldNames[iField].find("_y");
        if (found != string::npos) {
          writeFieldName = false;
          /*--- We have found a vector, so skip the Y component. ---*/
        }

        found = fieldNames[iField].find("_z");
        if (found != string::npos) {
          writeFieldName = false;
          /*--- We have found a vector, so skip the Z component. ---*/
        }

        if (writeFieldName) {
          string fieldDescription = "Results/Mesh/MultiStep_1/ResultDescription/" + thisFieldName;

          file.createDataSet<std::string>((fieldDescription + "/" + "DOFNames"), HighFive::DataSpace::From(DOFNames))
              .write(DOFNames);
          file.createDataSet<int>((fieldDescription + "/" + "DefinedOn"), HighFive::DataSpace(1)).write(1);
          file.createDataSet<std::string>((fieldDescription + "/" + "EntityNames"),
                                          HighFive::DataSpace::From(EntityName))
              .write(EntityName);
          file.createDataSet<int>((fieldDescription + "/" + "EntryType"), HighFive::DataSpace(1)).write(entryType);
          file.createDataSet<int>((fieldDescription + "/" + "NumDOFs"), HighFive::DataSpace(1)).write(entryType);
          file.createDataSet<int>((fieldDescription + "/" + "StepNumbers"), HighFive::DataSpace::From(stepNumbers))
              .write(stepNumbers);
          file.createDataSet<double>((fieldDescription + "/" + "StepValues"), HighFive::DataSpace::From(stepValues))
              .write(stepValues);
          file.createDataSet<std::string>((fieldDescription + "/" + "Unit"), HighFive::DataSpace::From(unit))
              .write(unit);
        }
      }
    }

    SU2_MPI::Barrier(SU2_MPI::GetComm());

    /*---Write the mesh in parallel ---*/

    HighFive::File file(H5_fileName, HighFive::File::ReadWrite | HighFive::File::OpenOrCreate,
                        HighFive::MPIOFileDriver(MPI_COMM_WORLD, MPI_INFO_NULL));
    HighFive::DataSet datasetMesh =
        file.createDataSet<double>(fullMeshPathName, HighFive::DataSpace(dataSorter->GetnPointsGlobal(), 3));
    datasetMesh.select({dataSorter->GetnPointCumulative(rank), 0}, {myPoint, 3}).write(twoDimBuffer);

    HighFive::DataSet dataSetRegionsNodesStatic =
        file.createDataSet<float>("Mesh/Nodes/Coordinates", HighFive::DataSpace(GlobalPoint, 3));
    dataSetRegionsNodesStatic.select({dataSorter->GetnPointCumulative(rank), 0}, {myPoint, 3}).write(twoDimBuffer);

    SU2_MPI::Barrier(SU2_MPI::GetComm());

    /*---Lambda to retrieve the element types and the offset to write the data from the current process ---*/
    auto copyToBuffer = [&](GEO_TYPE type, H5_Element_TYPE h5Type, unsigned long nElem, unsigned short nPoints) {
      for (iElem = 0; iElem < nElem; iElem++) {
        for (iNode = 0; iNode < nPoints; iNode++) {
          connBuf[iStorage + iNode] = int(dataSorter->GetElem_Connectivity(type, iElem, iNode) - 1);
          connectivity[iElemID][iNode] = int(dataSorter->GetElem_Connectivity(type, iElem, iNode));
        }

        iStorage += nPoints;
        offsetBuf[iElemID] = int(iStorage + dataSorter->GetnElemConnCumulative(rank));
        elementType[iElemID++] = h5Type;
      }
    };

    /*---Call the lambda for all possible cells ---*/
    copyToBuffer(LINE, H5LINE, nParallel_Line, N_POINTS_LINE);
    copyToBuffer(TRIANGLE, H5TRIANGLE, nParallel_Tria, N_POINTS_TRIANGLE);
    copyToBuffer(QUADRILATERAL, H5QUADRILATERAL, nParallel_Quad, N_POINTS_QUADRILATERAL);
    copyToBuffer(TETRAHEDRON, H5TETRAHEDRON, nParallel_Tetr, N_POINTS_TETRAHEDRON);
    copyToBuffer(HEXAHEDRON, H5HEXAHEDRON, nParallel_Hexa, N_POINTS_HEXAHEDRON);
    copyToBuffer(PRISM, H5PRISM, nParallel_Pris, N_POINTS_PRISM);
    copyToBuffer(PYRAMID, H5PYRAMID, nParallel_Pyra, N_POINTS_PYRAMID);

    string TopoPathName = "/Mesh/Elements/";
    string fullTopoStorPathName = TopoPathName + "TopoStorage";
    string fullTopoOffsetPathName = TopoPathName + "TopoOffset";
    string fullElementTypePathName = TopoPathName + "Types";

    HighFive::DataSet datasetTopoStr =
        file.createDataSet<int>(fullTopoStorPathName, HighFive::DataSpace(GlobalElemStorage));
    datasetTopoStr.select({dataSorter->GetnElemConnCumulative(rank)}, {myElemStorage}).write(connBuf);

    HighFive::DataSet datasetTopoOffSet =
        file.createDataSet<int>(fullTopoOffsetPathName, HighFive::DataSpace(GlobalElem));
    datasetTopoOffSet.select({dataSorter->GetnElemCumulative(rank)}, {myElem}).write(offsetBuf);

    HighFive::DataSet datasetElementType =
        file.createDataSet<int>(fullElementTypePathName, HighFive::DataSpace(GlobalElem));
    datasetElementType.select({dataSorter->GetnElemCumulative(rank)}, {myElem}).write(elementType);

    string fullTopoPathName = TopoPathName + "Connectivity";
    HighFive::DataSet datasetTopo =
        file.createDataSet<int>(fullTopoPathName, HighFive::DataSpace(GlobalElem, N_POINTS_HEXAHEDRON));

    datasetTopo.select({dataSorter->GetnElemCumulative(rank), 0}, {myElem, N_POINTS_HEXAHEDRON}).write(connectivity);
    SU2_MPI::Barrier(SU2_MPI::GetComm());
    setTopoState();
  }

  /*---for steady state simulations, replace the old solution with the latest one ---*/
  if (rank == MASTER_NODE && (getSteadyIteration() > 1)) {
    HighFive::File file(H5_fileName, HighFive::File::ReadWrite | HighFive::File::OpenOrCreate);

    string pathName = "/Results/Mesh/MultiStep_1/";
    string renamedStep = "Step_" + std::to_string(getSteadyIteration() - 1);
    string renamedPathName = pathName + renamedStep + "/";  //+thisFieldName+"/default/Nodes/Real";
    string steadyPathName = pathName + timestep + "/";      // thisFieldName+"/default/Nodes/Real";
    file.rename(steadyPathName, renamedPathName);
  }

  /*--- Loop over all variables that have been registered in the output and write the results. ---*/

  HighFive::File file(H5_fileName, HighFive::File::ReadWrite | HighFive::File::OpenOrCreate,
                      HighFive::MPIOFileDriver(MPI_COMM_WORLD, MPI_INFO_NULL));
  std::vector<std::vector<double> > vectorBuffer(myPoint, vector<double>(NCOORDS, 0.0));

  VarCounter = varStart;

  for (iField = varStart; iField < fieldNames.size(); iField++) {
    bool output_variable = true, isVector = false;
    size_t found = fieldNames[iField].find("_x");
    if (found != string::npos) {
      output_variable = true;
      isVector = true;
    }
    found = fieldNames[iField].find("_y");
    if (found != string::npos) {
      /*--- We have found a vector, so skip the Y component. ---*/
      output_variable = false;
      VarCounter++;
    }
    found = fieldNames[iField].find("_z");
    if (found != string::npos) {
      /*--- We have found a vector, so skip the Z component. ---*/
      output_variable = false;
      VarCounter++;
    }

    /*--- Write the point data as an <X,Y,Z> vector or a scalar. ---*/

    if (output_variable && isVector) {
      /*--- Load up the buffer for writing this rank's vector data. ---*/

      float val = 0.0;
      for (iPoint = 0; iPoint < myPoint; iPoint++) {
        for (iDim = 0; iDim < NCOORDS; iDim++) {
          if (nDim == 2 && iDim == 2) {
            dataBufferFloat[iPoint * NCOORDS + iDim] = 0.0;
            vectorBuffer[iPoint][iDim] = 0;

          } else {
            val = (float)dataSorter->GetData(VarCounter + iDim, iPoint);
            dataBufferFloat[iPoint * NCOORDS + iDim] = val;
            vectorBuffer[iPoint][iDim] = val;
          }
        }
      }

      string thisFieldName = fieldNames[iField];
      thisFieldName.erase(thisFieldName.end() - 2, thisFieldName.end());
      string pathName = "/Results/Mesh/MultiStep_1/";
      string fullPathName = pathName + timestep + "/" + thisFieldName + "/default/Nodes/Real";

      HighFive::DataSet dataset =
          file.createDataSet<double>(fullPathName, HighFive::DataSpace(dataSorter->GetnPointsGlobal(), NCOORDS));
      dataset.select({dataSorter->GetnPointCumulative(rank), 0}, {myPoint, NCOORDS}).write(vectorBuffer);

      VarCounter++;

    } else if (output_variable) {
      string dataSetName = fieldNames[iField];
      string pathName = "/Results/Mesh/MultiStep_1/";
      string fullPathName = pathName + timestep + "/" + dataSetName + "/default/Nodes/Real";
      HighFive::DataSet dataset =
          file.createDataSet<double>(fullPathName, HighFive::DataSpace(dataSorter->GetnPointsGlobal()));

      for (iPoint = 0; iPoint < myPoint; iPoint++) {
        float val = (float)dataSorter->GetData(VarCounter, iPoint);
        dataBufferFloat[iPoint] = val;
      }
      dataset.select({dataSorter->GetnPointCumulative(rank)}, {myPoint}).write(dataBufferFloat.data());
      VarCounter++;
    }
  }
}
