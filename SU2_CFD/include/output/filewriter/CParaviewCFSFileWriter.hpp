/*!
 * \file CParaviewXMLFileWriter.cpp
 * \brief Filewriter class for Paraview binary format.
 * \author T. Albring
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

#pragma once
#include "../../../../Common/include/CConfig.hpp"
#include "CFileWriter.hpp"

class CParaviewCFSFileWriter final : public CFileWriter {
 private:
  /*!
   * \brief Current timestep to be written
   */
  passivedouble timeStep;

  /*!
   * \brief The last timestep in an unsteady simulation to be written
   */
  unsigned long lastTimeStep;

  /*!
   * \brief Boolean storing whether mesh is already written to the H5 file
   */
  static bool topoState;

  /*!
   * \brief Int indicating how many intermediate steps are written to the disk
   */
  static int steadyIteration;

  /*!
   * \brief Boolean indicating if the simulation is steady-state
   */
  bool steadySim = NULL;

  /*!
   * \brief Int storing the current time iteration
   */
  int TimeIter;

 public:
  const static string fileExt;

  /*!
   * \brief Construct a file writer using field names and the data sorter.
   * \param[in] val_filename - The name of the file
   * \param[in] valDataSorter - The parallel sorted data to write
   * \param[in] config - The config file
   */
  CParaviewCFSFileWriter(string val_filename, CParallelDataSorter* valDataSorter, CConfig* config);

  /*!
   * \brief Construct a file writer using field names and the data sorter.
   * \param[in] valDataSorter - The parallel sorted data to write
   * \param[in] config - The config file
   */
  CParaviewCFSFileWriter(CParallelDataSorter* valDataSorter, CConfig* config);

  /*!
   * \brief Destructor
   */
  ~CParaviewCFSFileWriter() override;

  /*!
   * \brief Write sorted data to H5 file
   */
  void Write_Data(string val_filename) override;

  /*!
   * \brief Set the topoState true if the mesh is already written to the disk
   */
  int setTopoState() {
    topoState = true;
    return 0;
  };

  /*!
   * \brief Get the topoState to know whether or not the mesh is already written
   */
  bool getTopoState() { return topoState; };

  /*!
   * \brief Get the last timestep of the unsteady simulatiob
   */
  unsigned long getLastTimeStep() { return lastTimeStep; };

  /*!
   * \brief Get the current timestep
   */
  passivedouble getTime_Step() { return timeStep; };

  /*!
   * \brief determines if the simulation is steady state
   */
  bool isSteady() { return steadySim; };

  /*!
   * \brief increment the iteration by one for steady simulations
   */
  static void incrementSteadyIteration();

  /*!
   * \brief get the current number of intermediate steady solutions written to the disk
   */
  int getSteadyIteration() { return steadyIteration; };
};
