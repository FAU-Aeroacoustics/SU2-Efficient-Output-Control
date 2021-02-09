/*!
 * \file CFEM_DG_EulerSolver.cpp
 * \brief Main subroutines for solving finite element Euler flow problems
 * \author J. Alonso, E. van der Weide, T. Economon
 * \version 7.1.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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


#include "../../include/solvers/CFEM_DG_EulerSolver.hpp"
#include "../../../Common/include/toolboxes/printing_toolbox.hpp"
#include "../../include/fluid/CIdealGas.hpp"
#include "../../include/fluid/CVanDerWaalsGas.hpp"
#include "../../include/fluid/CPengRobinson.hpp"

CFEM_DG_EulerSolver::CFEM_DG_EulerSolver(void)
  : CFEM_DG_SolverBase() {}

CFEM_DG_EulerSolver::CFEM_DG_EulerSolver(CConfig        *config,
                                         unsigned short val_nDim,
                                         unsigned short iMesh)
  : CFEM_DG_SolverBase() {

  /*--- Store the multigrid level. ---*/
  MGLevel = iMesh;

  /*--- Dummy solver constructor that calls the SetNondim. routine in
        order to load the flow non-dim. information into the config class.
        This is needed to complete a partitioning for time-accurate local
        time stepping that depends on the flow state. ---*/
  nDim = val_nDim;
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  SetNondimensionalization(config, iMesh, false);
}

CFEM_DG_EulerSolver::CFEM_DG_EulerSolver(CGeometry      *geometry,
                                         CConfig        *config,
                                         unsigned short iMesh)
  : CFEM_DG_SolverBase(geometry, config, iMesh) {

  /*--- Set the gamma value and the number of variables. ---*/
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  nVar = nDim + 2;

  /*--- Define some auxiliary vectors related to the residual ---*/
  Residual_RMS = new su2double[nVar];     for(unsigned short iVar=0; iVar<nVar; ++iVar) Residual_RMS[iVar] = 1.e-35;
  Residual_Max = new su2double[nVar];     for(unsigned short iVar=0; iVar<nVar; ++iVar) Residual_Max[iVar] = 1.e-35;
  Point_Max    = new unsigned long[nVar]; for(unsigned short iVar=0; iVar<nVar; ++iVar) Point_Max[iVar]    = 0;

  Point_Max_Coord = new su2double*[nVar];
  for (unsigned short iVar=0; iVar<nVar; ++iVar) {
    Point_Max_Coord[iVar] = new su2double[nDim];
    for(unsigned short iDim=0; iDim<nDim; ++iDim) Point_Max_Coord[iVar][iDim] = 0.0;
  }

  /*--- Perform the non-dimensionalization for the flow equations using the
        specified reference values. ---*/
  SetNondimensionalization(config, iMesh, true);

  /*--- Check if we are executing a verification case. If so, the
        VerificationSolution object will be instantiated for a particular
        option from the available library of verification solutions. Note
        that this is done after SetNondim(), as problem-specific initial
        parameters are needed by the solution constructors. ---*/
  SetVerificationSolution(nDim, nVar, config);

  /*--- Read farfield conditions ---*/
  Density_Inf     = config->GetDensity_FreeStreamND();
  Pressure_Inf    = config->GetPressure_FreeStreamND();
  Velocity_Inf    = config->GetVelocity_FreeStreamND();
  Energy_Inf      = config->GetEnergy_FreeStreamND();
  Temperature_Inf = config->GetTemperature_FreeStreamND();
  Mach_Inf        = config->GetMach();

  /*--- Set the conservative variables of the free-stream. ---*/
  ConsVarFreeStream.resize(nVar);

  ConsVarFreeStream[0]      = Density_Inf;
  ConsVarFreeStream[nVar-1] = Density_Inf*Energy_Inf;

  for(unsigned short iDim=0; iDim<nDim; ++iDim)
    ConsVarFreeStream[iDim+1] = Density_Inf*Velocity_Inf[iDim];

  /*--- Set the entropy variables of the free-stream. ---*/
  EntropyVarFreeStream.resize(nVar);

  const su2double ovgm1 = 1.0/Gamma_Minus_One;
  const su2double s     = log(Pressure_Inf/pow(Density_Inf,Gamma));
  const su2double pInv  = 1.0/Pressure_Inf;

  su2double V2_Inf = 0.0;
  for(unsigned short iDim=0; iDim<nDim; ++iDim) {
    EntropyVarFreeStream[iDim+1] = Density_Inf*Velocity_Inf[iDim]*pInv;
    V2_Inf                      += Velocity_Inf[iDim]*Velocity_Inf[iDim];
  }

  EntropyVarFreeStream[0]      =  (Gamma-s)*ovgm1 - 0.5*pInv*Density_Inf*V2_Inf;
  EntropyVarFreeStream[nVar-1] = -pInv*Density_Inf;

  /*--- Parallel loop over the volume elements, if supperted. ---*/
  SU2_OMP_PARALLEL
  {
#ifdef HAVE_OMP
    const size_t omp_chunk_size_vol = computeStaticChunkSize(nVolElemTot, omp_get_num_threads(), 64);
#endif
    SU2_OMP_FOR_STAT(omp_chunk_size_vol)
    for(unsigned long i=0; i<nVolElemTot; ++i) {

      /*--- Allocate the memory to store the volume solution(s)
            and the residuals. ---*/
      volElem[i].AllocateCompressibleFlowVar(config, nVar);
      volElem[i].AllocateResiduals(config, nVar);

      /*--- Initialize the solution to a uniform flow. This is overruled
            when a restart takes place. ---*/
      volElem[i].SetConstantSolution(EntropyVarFreeStream.data(), nVar, 0);
    }
  }

  /*--- Determine the communication pattern. ---*/
  CMeshFEM_DG *DGGeometry = dynamic_cast<CMeshFEM_DG *>(geometry);
  if( !DGGeometry) SU2_MPI::Error(string("Dynamic cast failed"), CURRENT_FUNCTION);
  Prepare_MPI_Communication(DGGeometry, config);

  /*--- Set up the list of tasks to be carried out in the computational expensive
        part of the code. For the Runge-Kutta schemes this is typically the
        computation of the spatial residual, while for ADER this list contains
        the tasks to be done for one space time step. ---*/
  SetUpTaskList(config);
}

CFEM_DG_EulerSolver::~CFEM_DG_EulerSolver(void) {
}

void CFEM_DG_EulerSolver::SetNondimensionalization(CConfig        *config,
                                                   unsigned short iMesh,
                                                   const bool     writeOutput) {

  su2double Temperature_FreeStream   = 0.0, Mach2Vel_FreeStream   = 0.0, ModVel_FreeStream      = 0.0,
            Energy_FreeStream        = 0.0, ModVel_FreeStreamND   = 0.0, Velocity_Reynolds      = 0.0,
            Omega_FreeStream         = 0.0, Omega_FreeStreamND    = 0.0, Viscosity_FreeStream   = 0.0,
            Density_FreeStream       = 0.0, Pressure_FreeStream   = 0.0, Tke_FreeStream         = 0.0,
            Length_Ref               = 0.0, Density_Ref           = 0.0, Pressure_Ref           = 0.0,
            Velocity_Ref             = 0.0, Temperature_Ref       = 0.0, Time_Ref               = 0.0,
            Omega_Ref                = 0.0, Force_Ref             = 0.0, Gas_Constant_Ref       = 0.0,
            Viscosity_Ref            = 0.0, Conductivity_Ref      = 0.0, Energy_Ref             = 0.0,
            Froude                   = 0.0, Pressure_FreeStreamND = 0.0, Density_FreeStreamND   = 0.0,
            Temperature_FreeStreamND = 0.0, Gas_ConstantND        = 0.0, Viscosity_FreeStreamND = 0.0,
            Tke_FreeStreamND         = 0.0, Energy_FreeStreamND   = 0.0, Total_UnstTimeND       = 0.0,
            Delta_UnstTimeND         = 0.0;

  su2double Velocity_FreeStreamND[MAXNDIM] = {0.0};

  unsigned short iDim;

  /*--- Local variables ---*/
  su2double Alpha            = config->GetAoA()*PI_NUMBER/180.0;
  su2double Beta             = config->GetAoS()*PI_NUMBER/180.0;
  su2double Mach             = config->GetMach();
  su2double Reynolds         = config->GetReynolds();
  bool unsteady           = (config->GetTime_Marching() != NO);
  bool viscous            = config->GetViscous();
  bool grid_movement      = config->GetGrid_Movement();
  bool turbulent          = (config->GetKind_Solver() == FEM_RANS) || (config->GetKind_Solver() == FEM_LES);
  bool tkeNeeded          = ((turbulent) && ((config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST)));
  bool free_stream_temp   = (config->GetKind_FreeStreamOption() == TEMPERATURE_FS);
  bool reynolds_init      = (config->GetKind_InitOption() == REYNOLDS);

  /*--- Compute the Free Stream velocity, using the Mach number ---*/
  Pressure_FreeStream = config->GetPressure_FreeStream();
  Density_FreeStream  = config->GetDensity_FreeStream();
  Temperature_FreeStream  = config->GetTemperature_FreeStream();

  CFluidModel* auxFluidModel = nullptr;

  switch (config->GetKind_FluidModel()) {

    case STANDARD_AIR:

      if (config->GetSystemMeasurements() == SI) config->SetGas_Constant(287.058);
      else if (config->GetSystemMeasurements() == US) config->SetGas_Constant(1716.49);

      auxFluidModel = new CIdealGas(1.4, config->GetGas_Constant(), config->GetCompute_Entropy());
      if (free_stream_temp) {
        auxFluidModel->SetTDState_PT(Pressure_FreeStream, Temperature_FreeStream);
        Density_FreeStream = auxFluidModel->GetDensity();
        config->SetDensity_FreeStream(Density_FreeStream);
      }
      else {
        auxFluidModel->SetTDState_Prho(Pressure_FreeStream, Density_FreeStream );
        Temperature_FreeStream = auxFluidModel->GetTemperature();
        config->SetTemperature_FreeStream(Temperature_FreeStream);
      }
      break;

    case IDEAL_GAS:

      auxFluidModel = new CIdealGas(Gamma, config->GetGas_Constant(), config->GetCompute_Entropy());
      if (free_stream_temp) {
        auxFluidModel->SetTDState_PT(Pressure_FreeStream, Temperature_FreeStream);
        Density_FreeStream = auxFluidModel->GetDensity();
        config->SetDensity_FreeStream(Density_FreeStream);
      }
      else {
        auxFluidModel->SetTDState_Prho(Pressure_FreeStream, Density_FreeStream );
        Temperature_FreeStream = auxFluidModel->GetTemperature();
        config->SetTemperature_FreeStream(Temperature_FreeStream);
      }
      break;

    case VW_GAS:

      auxFluidModel = new CVanDerWaalsGas(Gamma, config->GetGas_Constant(),
                                          config->GetPressure_Critical(), config->GetTemperature_Critical());
      if (free_stream_temp) {
        auxFluidModel->SetTDState_PT(Pressure_FreeStream, Temperature_FreeStream);
        Density_FreeStream = auxFluidModel->GetDensity();
        config->SetDensity_FreeStream(Density_FreeStream);
      }
      else {
        auxFluidModel->SetTDState_Prho(Pressure_FreeStream, Density_FreeStream );
        Temperature_FreeStream = auxFluidModel->GetTemperature();
        config->SetTemperature_FreeStream(Temperature_FreeStream);
      }
      break;

    case PR_GAS:

      auxFluidModel = new CPengRobinson(Gamma, config->GetGas_Constant(), config->GetPressure_Critical(),
                                        config->GetTemperature_Critical(), config->GetAcentric_Factor());
      if (free_stream_temp) {
        auxFluidModel->SetTDState_PT(Pressure_FreeStream, Temperature_FreeStream);
        Density_FreeStream = auxFluidModel->GetDensity();
        config->SetDensity_FreeStream(Density_FreeStream);
      }
      else {
        auxFluidModel->SetTDState_Prho(Pressure_FreeStream, Density_FreeStream );
        Temperature_FreeStream = auxFluidModel->GetTemperature();
        config->SetTemperature_FreeStream(Temperature_FreeStream);
      }
      break;

  }

  Mach2Vel_FreeStream = auxFluidModel->GetSoundSpeed();

  /*--- Compute the Free Stream velocity, using the Mach number ---*/
  if (nDim == 2) {
    config->GetVelocity_FreeStream()[0] = cos(Alpha)*Mach*Mach2Vel_FreeStream;
    config->GetVelocity_FreeStream()[1] = sin(Alpha)*Mach*Mach2Vel_FreeStream;
  }
  if (nDim == 3) {
    config->GetVelocity_FreeStream()[0] = cos(Alpha)*cos(Beta)*Mach*Mach2Vel_FreeStream;
    config->GetVelocity_FreeStream()[1] = sin(Beta)*Mach*Mach2Vel_FreeStream;
    config->GetVelocity_FreeStream()[2] = sin(Alpha)*cos(Beta)*Mach*Mach2Vel_FreeStream;
  }

  /*--- Compute the modulus of the free stream velocity ---*/
  ModVel_FreeStream = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    ModVel_FreeStream += config->GetVelocity_FreeStream()[iDim]*config->GetVelocity_FreeStream()[iDim];
  ModVel_FreeStream = sqrt(ModVel_FreeStream); config->SetModVel_FreeStream(ModVel_FreeStream);

  /*--- Viscous initialization ---*/
  if (viscous) {

    /*--- The dimensional viscosity is needed to determine the free-stream conditions.
          To accomplish this, simply set the non-dimensional coefficients to the
          dimensional ones. This will be overruled later.---*/
    config->SetMu_RefND(config->GetMu_Ref());
    config->SetMu_Temperature_RefND(config->GetMu_Temperature_Ref());
    config->SetMu_SND(config->GetMu_S());

    config->SetMu_ConstantND(config->GetMu_Constant());

    /*--- Reynolds based initialization ---*/
    if (reynolds_init) {

      /*--- First, check if there is mesh motion. If yes, use the Mach
       number relative to the body to initialize the flow. ---*/
      if (grid_movement) Velocity_Reynolds = config->GetMach_Motion()*Mach2Vel_FreeStream;
      else Velocity_Reynolds = ModVel_FreeStream;

      /*--- For viscous flows, pressure will be computed from a density
            that is found from the Reynolds number. The viscosity is computed
            from the dimensional version of Sutherland's law or the constant
            viscosity, depending on the input option.---*/
      auxFluidModel->SetLaminarViscosityModel(config);

      Viscosity_FreeStream = auxFluidModel->GetLaminarViscosity();
      config->SetViscosity_FreeStream(Viscosity_FreeStream);

      Density_FreeStream = Reynolds*Viscosity_FreeStream/(Velocity_Reynolds*config->GetLength_Reynolds());
      config->SetDensity_FreeStream(Density_FreeStream);
      auxFluidModel->SetTDState_rhoT(Density_FreeStream, Temperature_FreeStream);
      Pressure_FreeStream = auxFluidModel->GetPressure();
      config->SetPressure_FreeStream(Pressure_FreeStream);
      Energy_FreeStream = auxFluidModel->GetStaticEnergy() + 0.5*ModVel_FreeStream*ModVel_FreeStream;

    }

    /*--- Thermodynamics quantities based initialization ---*/
    else {

      auxFluidModel->SetLaminarViscosityModel(config);
      Viscosity_FreeStream = auxFluidModel->GetLaminarViscosity();
      config->SetViscosity_FreeStream(Viscosity_FreeStream);
      Energy_FreeStream = auxFluidModel->GetStaticEnergy() + 0.5*ModVel_FreeStream*ModVel_FreeStream;

    }

    /*--- Turbulence kinetic energy ---*/
    Tke_FreeStream  = 3.0/2.0*(ModVel_FreeStream*ModVel_FreeStream*config->GetTurbulenceIntensity_FreeStream()*config->GetTurbulenceIntensity_FreeStream());

  }
  else {

    /*--- For inviscid flow, energy is calculated from the specified
     FreeStream quantities using the proper gas law. ---*/
    Energy_FreeStream = auxFluidModel->GetStaticEnergy() + 0.5*ModVel_FreeStream*ModVel_FreeStream;

  }

  /*-- Compute the freestream energy. ---*/
  if (tkeNeeded) { Energy_FreeStream += Tke_FreeStream; }; config->SetEnergy_FreeStream(Energy_FreeStream);

  /*--- Compute non dimensional quantities. By definition,
   Lref is one because we have converted the grid to meters. ---*/
  if (config->GetRef_NonDim() == DIMENSIONAL) {
    Pressure_Ref      = 1.0;
    Density_Ref       = 1.0;
    Temperature_Ref   = 1.0;
  }
  else if (config->GetRef_NonDim() == FREESTREAM_PRESS_EQ_ONE) {
    Pressure_Ref      = Pressure_FreeStream;     // Pressure_FreeStream = 1.0
    Density_Ref       = Density_FreeStream;      // Density_FreeStream = 1.0
    Temperature_Ref   = Temperature_FreeStream;  // Temperature_FreeStream = 1.0
  }
  else if (config->GetRef_NonDim() == FREESTREAM_VEL_EQ_MACH) {
    Pressure_Ref      = Gamma*Pressure_FreeStream; // Pressure_FreeStream = 1.0/Gamma
    Density_Ref       = Density_FreeStream;        // Density_FreeStream = 1.0
    Temperature_Ref   = Temperature_FreeStream;    // Temp_FreeStream = 1.0
  }
  else if (config->GetRef_NonDim() == FREESTREAM_VEL_EQ_ONE) {
    Pressure_Ref      = Mach*Mach*Gamma*Pressure_FreeStream; // Pressure_FreeStream = 1.0/(Gamma*(M_inf)^2)
    Density_Ref       = Density_FreeStream;        // Density_FreeStream = 1.0
    Temperature_Ref   = Temperature_FreeStream;    // Temp_FreeStream = 1.0
  }
  config->SetPressure_Ref(Pressure_Ref);
  config->SetDensity_Ref(Density_Ref);
  config->SetTemperature_Ref(Temperature_Ref);

  Length_Ref        = 1.0;                                                         config->SetLength_Ref(Length_Ref);
  Velocity_Ref      = sqrt(config->GetPressure_Ref()/config->GetDensity_Ref());    config->SetVelocity_Ref(Velocity_Ref);
  Time_Ref          = Length_Ref/Velocity_Ref;                                     config->SetTime_Ref(Time_Ref);
  Omega_Ref         = Velocity_Ref/Length_Ref;                                     config->SetOmega_Ref(Omega_Ref);
  Force_Ref         = Velocity_Ref*Velocity_Ref/Length_Ref;                        config->SetForce_Ref(Force_Ref);
  Gas_Constant_Ref  = Velocity_Ref*Velocity_Ref/config->GetTemperature_Ref();      config->SetGas_Constant_Ref(Gas_Constant_Ref);
  Viscosity_Ref     = config->GetDensity_Ref()*Velocity_Ref*Length_Ref;            config->SetViscosity_Ref(Viscosity_Ref);
  Conductivity_Ref  = Viscosity_Ref*Gas_Constant_Ref;                              config->SetConductivity_Ref(Conductivity_Ref);
  Froude            = ModVel_FreeStream/sqrt(STANDARD_GRAVITY*Length_Ref);         config->SetFroude(Froude);

  /*--- Divide by reference values, to compute the non-dimensional free-stream values ---*/
  Pressure_FreeStreamND = Pressure_FreeStream/config->GetPressure_Ref();  config->SetPressure_FreeStreamND(Pressure_FreeStreamND);
  Density_FreeStreamND  = Density_FreeStream/config->GetDensity_Ref();    config->SetDensity_FreeStreamND(Density_FreeStreamND);

  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_FreeStreamND[iDim] = config->GetVelocity_FreeStream()[iDim]/Velocity_Ref; config->SetVelocity_FreeStreamND(Velocity_FreeStreamND[iDim], iDim);
  }

  Temperature_FreeStreamND = Temperature_FreeStream/config->GetTemperature_Ref(); config->SetTemperature_FreeStreamND(Temperature_FreeStreamND);

  Gas_ConstantND = config->GetGas_Constant()/Gas_Constant_Ref;    config->SetGas_ConstantND(Gas_ConstantND);


  ModVel_FreeStreamND = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) ModVel_FreeStreamND += Velocity_FreeStreamND[iDim]*Velocity_FreeStreamND[iDim];
  ModVel_FreeStreamND    = sqrt(ModVel_FreeStreamND); config->SetModVel_FreeStreamND(ModVel_FreeStreamND);

  Viscosity_FreeStreamND = Viscosity_FreeStream / Viscosity_Ref;   config->SetViscosity_FreeStreamND(Viscosity_FreeStreamND);

  Tke_FreeStream  = 3.0/2.0*(ModVel_FreeStream*ModVel_FreeStream*config->GetTurbulenceIntensity_FreeStream()*config->GetTurbulenceIntensity_FreeStream());
  config->SetTke_FreeStream(Tke_FreeStream);

  Tke_FreeStreamND  = 3.0/2.0*(ModVel_FreeStreamND*ModVel_FreeStreamND*config->GetTurbulenceIntensity_FreeStream()*config->GetTurbulenceIntensity_FreeStream());
  config->SetTke_FreeStreamND(Tke_FreeStreamND);

  Omega_FreeStream = Density_FreeStream*Tke_FreeStream/max((Viscosity_FreeStream*config->GetTurb2LamViscRatio_FreeStream()), 1.e-25);
  config->SetOmega_FreeStream(Omega_FreeStream);

  Omega_FreeStreamND = Density_FreeStreamND*Tke_FreeStreamND/max((Viscosity_FreeStreamND*config->GetTurb2LamViscRatio_FreeStream()), 1.e-25);
  config->SetOmega_FreeStreamND(Omega_FreeStreamND);

  /*--- Set viscosity ND constants before defining the visc. model of the fluid models. ---*/

  if (viscous) {
    /*--- Constant viscosity model. ---*/
    config->SetMu_ConstantND(config->GetMu_Constant()/Viscosity_Ref);

    /*--- Sutherland's model. ---*/
    config->SetMu_RefND(config->GetMu_Ref()/Viscosity_Ref);
    config->SetMu_SND(config->GetMu_S()/config->GetTemperature_Ref());
    config->SetMu_Temperature_RefND(config->GetMu_Temperature_Ref()/config->GetTemperature_Ref());

    /*--- Constant thermal conductivity model. ---*/
    config->SetKt_ConstantND(config->GetKt_Constant()/Conductivity_Ref);
  }
  

  /*--- Initialize the dimensionless Fluid Model that will be used to solve the dimensionless problem ---*/

  /*--- Auxilary (dimensional) FluidModel no longer needed. ---*/
  delete auxFluidModel;

  /*--- Create one final fluid model object per OpenMP thread to be able to use them in parallel.
   *    GetFluidModel() should be used to automatically access the "right" object of each thread. ---*/

  assert(FluidModel.empty() && "Potential memory leak!");
  FluidModel.resize(omp_get_max_threads());

  SU2_OMP_PARALLEL
  {
    const int thread = omp_get_thread_num();

    switch (config->GetKind_FluidModel()) {

      case STANDARD_AIR:
        FluidModel[thread] = new CIdealGas(1.4, Gas_ConstantND, config->GetCompute_Entropy());
        break;

      case IDEAL_GAS:
        FluidModel[thread] = new CIdealGas(Gamma, Gas_ConstantND, config->GetCompute_Entropy());
        break;

      case VW_GAS:
        FluidModel[thread] = new CVanDerWaalsGas(Gamma, Gas_ConstantND,
                                                 config->GetPressure_Critical() /config->GetPressure_Ref(),
                                                 config->GetTemperature_Critical()/config->GetTemperature_Ref());
        break;

      case PR_GAS:
        FluidModel[thread] = new CPengRobinson(Gamma, Gas_ConstantND,
                                               config->GetPressure_Critical() /config->GetPressure_Ref(),
                                               config->GetTemperature_Critical()/config->GetTemperature_Ref(),
                                               config->GetAcentric_Factor());
        break;
    }

    GetFluidModel()->SetEnergy_Prho(Pressure_FreeStreamND, Density_FreeStreamND);
    if (viscous) {
      GetFluidModel()->SetLaminarViscosityModel(config);
      GetFluidModel()->SetThermalConductivityModel(config);
    }

  } // end SU2_OMP_PARALLEL

  Energy_FreeStreamND = GetFluidModel()->GetStaticEnergy() + 0.5*ModVel_FreeStreamND*ModVel_FreeStreamND;

  if (tkeNeeded) { Energy_FreeStreamND += Tke_FreeStreamND; };  config->SetEnergy_FreeStreamND(Energy_FreeStreamND);

  Energy_Ref = Energy_FreeStream/Energy_FreeStreamND; config->SetEnergy_Ref(Energy_Ref);

  Total_UnstTimeND = config->GetTotal_UnstTime() / Time_Ref;    config->SetTotal_UnstTimeND(Total_UnstTimeND);
  Delta_UnstTimeND = config->GetDelta_UnstTime() / Time_Ref;    config->SetDelta_UnstTimeND(Delta_UnstTimeND);

  /*--- Write output to the console if this is required and if this is the master node and first domain ---*/
  if ((rank == MASTER_NODE) && (iMesh == MESH_0) && writeOutput) {

    cout.precision(6);

    if (viscous) {
      cout << "Viscous flow: Computing pressure using the ideal gas law" << endl;
      cout << "based on the free-stream temperature and a density computed" << endl;
      cout << "from the Reynolds number." << endl;
    } else {
      cout << "Inviscid flow: Computing density based on free-stream" << endl;
      cout << "temperature and pressure using the ideal gas law." << endl;
    }

    if (grid_movement) cout << "Force coefficients computed using MACH_MOTION." << endl;
    else cout << "Force coefficients computed using free-stream values." << endl;

    stringstream NonDimTableOut, ModelTableOut;
    stringstream Unit;

    cout << endl;
    PrintingToolbox::CTablePrinter ModelTable(&ModelTableOut);
    ModelTableOut <<"-- Models:"<< endl;

    ModelTable.AddColumn("Viscosity Model", 25);
    ModelTable.AddColumn("Conductivity Model", 26);
    ModelTable.AddColumn("Fluid Model", 25);
    ModelTable.SetAlign(PrintingToolbox::CTablePrinter::RIGHT);
    ModelTable.PrintHeader();

    PrintingToolbox::CTablePrinter NonDimTable(&NonDimTableOut);
    NonDimTable.AddColumn("Name", 22);
    NonDimTable.AddColumn("Dim. value", 14);
    NonDimTable.AddColumn("Ref. value", 14);
    NonDimTable.AddColumn("Unit", 10);
    NonDimTable.AddColumn("Non-dim. value", 14);
    NonDimTable.SetAlign(PrintingToolbox::CTablePrinter::RIGHT);

    NonDimTableOut <<"-- Fluid properties:"<< endl;

    NonDimTable.PrintHeader();

    if (viscous) {

      switch(config->GetKind_ViscosityModel()){
      case CONSTANT_VISCOSITY:
        ModelTable << "CONSTANT_VISCOSITY";
        if      (config->GetSystemMeasurements() == SI) Unit << "N.s/m^2";
        else if (config->GetSystemMeasurements() == US) Unit << "lbf.s/ft^2";
        NonDimTable << "Viscosity" << config->GetMu_Constant() << config->GetMu_Constant()/config->GetMu_ConstantND() << Unit.str() << config->GetMu_ConstantND();
        Unit.str("");
        NonDimTable.PrintFooter();
        break;

      case SUTHERLAND:
        ModelTable << "SUTHERLAND";
        if      (config->GetSystemMeasurements() == SI) Unit << "N.s/m^2";
        else if (config->GetSystemMeasurements() == US) Unit << "lbf.s/ft^2";
        NonDimTable << "Ref. Viscosity" <<  config->GetMu_Ref() <<  config->GetViscosity_Ref() << Unit.str() << config->GetMu_RefND();
        Unit.str("");
        if      (config->GetSystemMeasurements() == SI) Unit << "K";
        else if (config->GetSystemMeasurements() == US) Unit << "R";
        NonDimTable << "Sutherland Temp." << config->GetMu_Temperature_Ref() <<  config->GetTemperature_Ref() << Unit.str() << config->GetMu_Temperature_RefND();
        Unit.str("");
        if      (config->GetSystemMeasurements() == SI) Unit << "K";
        else if (config->GetSystemMeasurements() == US) Unit << "R";
        NonDimTable << "Sutherland Const." << config->GetMu_S() << config->GetTemperature_Ref() << Unit.str() << config->GetMu_SND();
        Unit.str("");
        NonDimTable.PrintFooter();
        break;

      }
      switch(config->GetKind_ConductivityModel()){
      case CONSTANT_PRANDTL:
        ModelTable << "CONSTANT_PRANDTL";
        NonDimTable << "Prandtl (Lam.)"  << "-" << "-" << "-" << config->GetPrandtl_Lam();
        Unit.str("");
        NonDimTable << "Prandtl (Turb.)" << "-" << "-" << "-" << config->GetPrandtl_Turb();
        Unit.str("");
        NonDimTable.PrintFooter();
        break;

      case CONSTANT_CONDUCTIVITY:
        ModelTable << "CONSTANT_CONDUCTIVITY";
        Unit << "W/m^2.K";
        NonDimTable << "Molecular Cond." << config->GetKt_Constant() << config->GetKt_Constant()/config->GetKt_ConstantND() << Unit.str() << config->GetKt_ConstantND();
        Unit.str("");
        NonDimTable.PrintFooter();
        break;

      }
    } else {
      ModelTable << "-" << "-";
    }

    if      (config->GetSystemMeasurements() == SI) Unit << "N.m/kg.K";
    else if (config->GetSystemMeasurements() == US) Unit << "lbf.ft/slug.R";
    NonDimTable << "Gas Constant" << config->GetGas_Constant() << config->GetGas_Constant_Ref() << Unit.str() << config->GetGas_ConstantND();
    Unit.str("");
    if      (config->GetSystemMeasurements() == SI) Unit << "N.m/kg.K";
    else if (config->GetSystemMeasurements() == US) Unit << "lbf.ft/slug.R";
    NonDimTable << "Spec. Heat Ratio" << "-" << "-" << "-" << Gamma;
    Unit.str("");

    switch(config->GetKind_FluidModel()){
    case STANDARD_AIR:
      ModelTable << "STANDARD_AIR";
      break;
    case IDEAL_GAS:
      ModelTable << "IDEAL_GAS";
      break;
    case VW_GAS:
      ModelTable << "VW_GAS";
      break;
    case PR_GAS:
      ModelTable << "PR_GAS";
      break;
    }

    if (config->GetKind_FluidModel() == VW_GAS || config->GetKind_FluidModel() == PR_GAS){
        NonDimTable << "Critical Pressure" << config->GetPressure_Critical() << config->GetPressure_Ref() << Unit.str() << config->GetPressure_Critical() /config->GetPressure_Ref();
        Unit.str("");
        Unit << "K";
        NonDimTable << "Critical Temperature" << config->GetTemperature_Critical() << config->GetTemperature_Ref() << Unit.str() << config->GetTemperature_Critical() /config->GetTemperature_Ref();
        Unit.str("");
    }
    NonDimTable.PrintFooter();

    NonDimTableOut <<"-- Initial and free-stream conditions:"<< endl;

    NonDimTable.PrintHeader();

    if      (config->GetSystemMeasurements() == SI) Unit << "Pa";
    else if (config->GetSystemMeasurements() == US) Unit << "psf";
    NonDimTable << "Static Pressure" << config->GetPressure_FreeStream() << config->GetPressure_Ref() << Unit.str() << config->GetPressure_FreeStreamND();
    Unit.str("");
    if      (config->GetSystemMeasurements() == SI) Unit << "kg/m^3";
    else if (config->GetSystemMeasurements() == US) Unit << "slug/ft^3";
    NonDimTable << "Density" << config->GetDensity_FreeStream() << config->GetDensity_Ref() << Unit.str() << config->GetDensity_FreeStreamND();
    Unit.str("");
    if      (config->GetSystemMeasurements() == SI) Unit << "K";
    else if (config->GetSystemMeasurements() == US) Unit << "R";
    NonDimTable << "Temperature" << config->GetTemperature_FreeStream() << config->GetTemperature_Ref() << Unit.str() << config->GetTemperature_FreeStreamND();
    Unit.str("");
    if      (config->GetSystemMeasurements() == SI) Unit << "m^2/s^2";
    else if (config->GetSystemMeasurements() == US) Unit << "ft^2/s^2";
    NonDimTable << "Total Energy" << config->GetEnergy_FreeStream() << config->GetEnergy_Ref() << Unit.str() << config->GetEnergy_FreeStreamND();
    Unit.str("");
    if      (config->GetSystemMeasurements() == SI) Unit << "m/s";
    else if (config->GetSystemMeasurements() == US) Unit << "ft/s";
    NonDimTable << "Velocity-X" << config->GetVelocity_FreeStream()[0] << config->GetVelocity_Ref() << Unit.str() << config->GetVelocity_FreeStreamND()[0];
    NonDimTable << "Velocity-Y" << config->GetVelocity_FreeStream()[1] << config->GetVelocity_Ref() << Unit.str() << config->GetVelocity_FreeStreamND()[1];
    if (nDim == 3) {
      NonDimTable << "Velocity-Z" << config->GetVelocity_FreeStream()[2] << config->GetVelocity_Ref() << Unit.str() << config->GetVelocity_FreeStreamND()[2];
    }
    NonDimTable << "Velocity Magnitude" << config->GetModVel_FreeStream() << config->GetVelocity_Ref() << Unit.str() << config->GetModVel_FreeStreamND();
    Unit.str("");

    if (viscous) {
      NonDimTable.PrintFooter();
      if      (config->GetSystemMeasurements() == SI) Unit << "N.s/m^2";
      else if (config->GetSystemMeasurements() == US) Unit << "lbf.s/ft^2";
      NonDimTable << "Viscosity" << config->GetViscosity_FreeStream() << config->GetViscosity_Ref() << Unit.str() << config->GetViscosity_FreeStreamND();
      Unit.str("");
      if      (config->GetSystemMeasurements() == SI) Unit << "W/m^2.K";
      else if (config->GetSystemMeasurements() == US) Unit << "lbf/ft.s.R";
      NonDimTable << "Conductivity" << "-" << config->GetConductivity_Ref() << Unit.str() << "-";
      Unit.str("");
      if (turbulent) {
        if      (config->GetSystemMeasurements() == SI) Unit << "m^2/s^2";
        else if (config->GetSystemMeasurements() == US) Unit << "ft^2/s^2";
        NonDimTable << "Turb. Kin. Energy" << config->GetTke_FreeStream() << config->GetTke_FreeStream()/config->GetTke_FreeStreamND() << Unit.str() << config->GetTke_FreeStreamND();
        Unit.str("");
        if      (config->GetSystemMeasurements() == SI) Unit << "1/s";
        else if (config->GetSystemMeasurements() == US) Unit << "1/s";
        NonDimTable << "Spec. Dissipation" << config->GetOmega_FreeStream() << config->GetOmega_FreeStream()/config->GetOmega_FreeStreamND() << Unit.str() << config->GetOmega_FreeStreamND();
        Unit.str("");
      }
    }

    NonDimTable.PrintFooter();
    NonDimTable << "Mach Number" << "-" << "-" << "-" << config->GetMach();
    if (viscous) {
      NonDimTable << "Reynolds Number" << "-" << "-" << "-" << config->GetReynolds();
    }
    NonDimTable.PrintFooter();
    ModelTable.PrintFooter();

    if (unsteady){
      NonDimTableOut << "-- Unsteady conditions" << endl;
      NonDimTable.PrintHeader();
      NonDimTable << "Total Time" << config->GetMax_Time() << config->GetTime_Ref() << "s" << config->GetMax_Time()/config->GetTime_Ref();
      Unit.str("");
      NonDimTable << "Time Step" << config->GetTime_Step() << config->GetTime_Ref() << "s" << config->GetDelta_UnstTimeND();
      Unit.str("");
      NonDimTable.PrintFooter();
    }

    cout << ModelTableOut.str();
    cout << NonDimTableOut.str();
  }
}

void CFEM_DG_EulerSolver::SetUpTaskList(CConfig *config) {

  /* Check whether an ADER space-time step must be carried out.
     When only a spatial Jacobian is computed this is false per definition.  */
  if( (config->GetKind_TimeIntScheme_Flow() == ADER_DG) &&
     !(config->GetJacobian_Spatial_Discretization_Only()) ) {

    /*------------------------------------------------------------------------*/
    /* ADER time integration with local time stepping. There are 2^(M-1)      */
    /* subtime steps, where M is the number of time levels. For each of       */
    /* the subtime steps a number of tasks must be carried out, hence the     */
    /* total list can become rather lengthy.                                  */
    /*------------------------------------------------------------------------*/

    /*--- Determine whether or not the predictor solution must be interpolated
          for the time levels on this rank. Make a distinction between owned
          and halo elements. Also determine whether or not boundary conditions
          are present and whether or not these depend on halo data. ---*/
    const unsigned short nTimeLevels = config->GetnLevels_TimeAccurateLTS();
    vector<bool> interpolOwnedElem(nTimeLevels, false);
    vector<bool> interpolHaloElem(nTimeLevels, false);
    vector<bool> BCPresent(nTimeLevels, false);
    vector<bool> BCDependOnHalos(nTimeLevels, false);

    for(unsigned short level=0; level<nTimeLevels; ++level) {

      /* First check the boundary conditions. */
      for(unsigned short iMarker=0; iMarker<config->GetnMarker_All(); iMarker++) {

        const unsigned long surfElemBeg = boundaries[iMarker].nSurfElem[level];
        const unsigned long surfElemEnd = boundaries[iMarker].nSurfElem[level+1];
        if(surfElemEnd > surfElemBeg) {
          BCPresent[level] = true;
          if( boundaries[iMarker].haloInfoNeededForBC ) BCDependOnHalos[level] = true;
        }
      }

      /* Determine whether or not the owned elements must be interpolated. */
      if(nVolElemOwnedPerTimeLevel[level+1] > nVolElemOwnedPerTimeLevel[level])
        interpolOwnedElem[level] = true;
      if(nMatchingInternalFacesLocalElem[level+1] > nMatchingInternalFacesLocalElem[level])
       interpolOwnedElem[level] = true;
      if(nMatchingInternalFacesWithHaloElem[level+1] > nMatchingInternalFacesWithHaloElem[level])
        interpolOwnedElem[level] = true;
      if( BCPresent[level] )
        interpolOwnedElem[level] = true;

      /* Determine whether or not the halo elements must be interpolated. */
      if(nMatchingInternalFacesWithHaloElem[level+1] > nMatchingInternalFacesWithHaloElem[level])
        interpolHaloElem[level] = true;
      if( BCDependOnHalos[level] )
        interpolHaloElem[level] = true;
    }

    /* Determine the number of subtime steps and abbreviate the number of
       time integration points a bit easier.  */
    const unsigned int   nSubTimeSteps          = pow(2,nTimeLevels-1);
    const unsigned short nTimeIntegrationPoints = config->GetnTimeIntegrationADER_DG();

    /* Define the two dimensional vector to store the latest index for a
       certain task for every time level. */
    vector<vector<int> > indexInList(CTaskDefinition::ADER_UPDATE_SOLUTION+1,
                                     vector<int>(nTimeLevels, -1));

    /* Loop over the number of subtime steps in the algorithm. */
    for(unsigned int step=0; step<nSubTimeSteps; ++step) {

      /* Determine the time level for which an update must be carried out
         for this step. */
      unsigned int ii = step+1;
      unsigned short timeLevel = 0;
      while( !(ii%2) ) {ii/=2; ++timeLevel;}

      /* Definition of the variable to store previous indices of
         tasks that must have been completed. */
      int prevInd[5];

      /* Carry out the predictor step of the communication elements of level 0
         if these elements are present on this rank. */
      unsigned long elemBeg = nVolElemOwnedPerTimeLevel[0]
                            + nVolElemInternalPerTimeLevel[0];
      unsigned long elemEnd = nVolElemOwnedPerTimeLevel[1];
      if(elemEnd > elemBeg) {
        prevInd[0] = indexInList[CTaskDefinition::ADER_UPDATE_SOLUTION][0];
        indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_COMM_ELEMENTS][0] = (int)tasksList.size();
        tasksList.push_back(CTaskDefinition(CTaskDefinition::ADER_PREDICTOR_STEP_COMM_ELEMENTS, 0, prevInd[0]));
      }

      /* Initiate the communication of elements of level 0, if there is
         something to be communicated. */
#ifdef HAVE_MPI
      if( commRequests[0].size() ) {
        if(elemEnd > elemBeg)   // Data needs to be computed before sending.
          prevInd[0] = indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_COMM_ELEMENTS][0];
        else                    // Data only needs to be received.
          prevInd[0] = indexInList[CTaskDefinition::ADER_TIME_INTERPOLATE_HALO_ELEMENTS][0];

        /* Make sure that any previous communication has been completed. */
        prevInd[1] = indexInList[CTaskDefinition::COMPLETE_REVERSE_MPI_COMMUNICATION][0];

        /* Create the task. */
        indexInList[CTaskDefinition::INITIATE_MPI_COMMUNICATION][0] = tasksList.size();
        tasksList.push_back(CTaskDefinition(CTaskDefinition::INITIATE_MPI_COMMUNICATION, 0,
                                            prevInd[0], prevInd[1]));
      }
#endif

      /* Check if the time level is not nTimeLevels-1. In that case carry out
         the predictor step of the communication elements for the next time
         level and initiate their communication, if appropriate on this rank. */
      if(timeLevel < (nTimeLevels-1)) {
        const unsigned short nL = timeLevel+1;

        /* Carry out the predictor step of the communication elements of level nL
           if these elements are present on this rank. */
        elemBeg = nVolElemOwnedPerTimeLevel[nL] + nVolElemInternalPerTimeLevel[nL];
        elemEnd = nVolElemOwnedPerTimeLevel[nL+1];

        if(elemEnd > elemBeg) {
          prevInd[0] = indexInList[CTaskDefinition::ADER_UPDATE_SOLUTION][nL];
          indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_COMM_ELEMENTS][nL] = (int)tasksList.size();
          tasksList.push_back(CTaskDefinition(CTaskDefinition::ADER_PREDICTOR_STEP_COMM_ELEMENTS, nL, prevInd[0]));
        }

        /* Initiate the communication of elements of level nL, if there is
           something to be communicated. */
#ifdef HAVE_MPI
        if( commRequests[nL].size() ) {
          if(elemEnd > elemBeg)   // Data needs to be computed before sending.
            prevInd[0] = indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_COMM_ELEMENTS][nL];
          else                    // Data only needs to be received.
            prevInd[0] = indexInList[CTaskDefinition::ADER_TIME_INTERPOLATE_HALO_ELEMENTS][nL];

          /* Make sure that any previous communication has been completed. */
          prevInd[1] = indexInList[CTaskDefinition::COMPLETE_REVERSE_MPI_COMMUNICATION][nL];

          /* Create the actual task. */
          indexInList[CTaskDefinition::INITIATE_MPI_COMMUNICATION][nL] = tasksList.size();
          tasksList.push_back(CTaskDefinition(CTaskDefinition::INITIATE_MPI_COMMUNICATION, nL,
                                              prevInd[0], prevInd[1]));
        }
#endif
      }

      /* Carry out the predictor step of the internal elements of time level 0,
         if these elements are present on this rank. */
      elemBeg = nVolElemOwnedPerTimeLevel[0];
      elemEnd = nVolElemOwnedPerTimeLevel[0] + nVolElemInternalPerTimeLevel[0];

      if(elemEnd > elemBeg) {
        prevInd[0] = indexInList[CTaskDefinition::ADER_UPDATE_SOLUTION][0];
        indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_INTERNAL_ELEMENTS][0] = (int)tasksList.size();
        tasksList.push_back(CTaskDefinition(CTaskDefinition::ADER_PREDICTOR_STEP_INTERNAL_ELEMENTS, 0, prevInd[0]));
      }

      /* Determine the tasks to be completed before the communication of time
         level 0 can be completed. */
      prevInd[0] = prevInd[1] = prevInd[2] = -1;
#ifdef HAVE_MPI
      if( commRequests[0].size() )
        prevInd[0] = indexInList[CTaskDefinition::INITIATE_MPI_COMMUNICATION][0];
#endif

      if( elementsSendSelfComm[0].size() ) {
        prevInd[1] = indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_COMM_ELEMENTS][0];
        prevInd[2] = indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_INTERNAL_ELEMENTS][0];
      }

      /* Make sure that the -1 in prevInd, if any, are numbered last. */
      sort(prevInd, prevInd+3, greater<int>());

      /* Complete the communication of time level 0, if there is something
         to be completed. */
      if(prevInd[0] > -1) {
        indexInList[CTaskDefinition::COMPLETE_MPI_COMMUNICATION][0] = (int)tasksList.size();
        tasksList.push_back(CTaskDefinition(CTaskDefinition::COMPLETE_MPI_COMMUNICATION, 0,
                                            prevInd[0], prevInd[1], prevInd[2]));
      }

      /* Check if the time level is not nTimeLevels-1. In that case carry out
         the predictor step of the internal elements for the next time
         level and complete the communication for this time level,
         if appropriate on this rank. */
      if(timeLevel < (nTimeLevels-1)) {
        const unsigned short nL = timeLevel+1;

        /* Carry out the predictor step of the internal elements of level nL
           if these elements are present on this rank. */
        elemBeg = nVolElemOwnedPerTimeLevel[nL];
        elemEnd = nVolElemOwnedPerTimeLevel[nL] + nVolElemInternalPerTimeLevel[nL];

        if(elemEnd > elemBeg) {
          prevInd[0] = indexInList[CTaskDefinition::ADER_UPDATE_SOLUTION][nL];
          indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_INTERNAL_ELEMENTS][nL] = (int)tasksList.size();
          tasksList.push_back(CTaskDefinition(CTaskDefinition::ADER_PREDICTOR_STEP_INTERNAL_ELEMENTS, nL, prevInd[0]));
        }

        /* Determine the tasks to be completed before the communication of time
           level nL can be completed. */
        prevInd[0] = prevInd[1] = prevInd[2] = -1;
#ifdef HAVE_MPI
        if( commRequests[nL].size() )
          prevInd[0] = indexInList[CTaskDefinition::INITIATE_MPI_COMMUNICATION][nL];
#endif

        if( elementsSendSelfComm[nL].size() ) {
          prevInd[1] = indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_COMM_ELEMENTS][nL];
          prevInd[2] = indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_INTERNAL_ELEMENTS][nL];
        }

        /* Make sure that the -1 in prevInd, if any, are numbered last. */
        sort(prevInd, prevInd+3, greater<int>());

        /* Complete the communication of time level nL, if there is something
           to be completed. */
        if(prevInd[0] > -1) {
          indexInList[CTaskDefinition::COMPLETE_MPI_COMMUNICATION][nL] = (int)tasksList.size();
          tasksList.push_back(CTaskDefinition(CTaskDefinition::COMPLETE_MPI_COMMUNICATION, nL,
                                              prevInd[0], prevInd[1], prevInd[2]));
        }
      }

      /* Loop over the time integration points to compute the space time
         integral for the corrector step. */
      for(unsigned short intPoint=0; intPoint<nTimeIntegrationPoints; ++intPoint) {

        /* Loop over the time levels to be treated. */
        for(unsigned short level=0; level<=timeLevel; ++level) {

          /* Easier storage of the number of owned elements for this level
             and the number of owned elements of the adjacent time level
             that are needed for the computation of the surface integral. */
          unsigned long nOwnedElem = nVolElemOwnedPerTimeLevel[level+1]
                                   - nVolElemOwnedPerTimeLevel[level];

          unsigned long nAdjOwnedElem = 0;
          if(level < (nTimeLevels-1))
            nAdjOwnedElem = ownedElemAdjLowTimeLevel[level+1].size();

          /* Check if the solution of the owned elements of the current time
             level and the adjacent owned elements of the next time level must
             be interpolated to the integration point intPoint. */
          if( interpolOwnedElem[level] ) {

            /* Determine the tasks that should have been completed before this
               task can be carried out. This depends on the time integration
               point. */
            if(intPoint == 0) {

              /* First time integration point. The predictor steps of the
                 owned elements must have been completed before this task
                 can be carried out. */
              if( nOwnedElem ) {
                prevInd[0] = indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_COMM_ELEMENTS][level];
                prevInd[1] = indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_INTERNAL_ELEMENTS][level];
              }
              else {
                prevInd[0] = prevInd[1] = -1;
              }

              if( nAdjOwnedElem ) {
                prevInd[2] = indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_COMM_ELEMENTS][level+1];
                prevInd[3] = indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_INTERNAL_ELEMENTS][level+1];
              }
              else {
                prevInd[2] = prevInd[3] = -1;
              }

              prevInd[4] = -1;
            }
            else {

              /* Not the first integration point. The residual computations of
                 the previous integration point must have been completed, because
                 the solution in the work vectors is overwritten. */
              prevInd[0] = indexInList[CTaskDefinition::VOLUME_RESIDUAL][level];
              prevInd[1] = indexInList[CTaskDefinition::BOUNDARY_CONDITIONS_DEPEND_ON_OWNED][level];
              prevInd[2] = indexInList[CTaskDefinition::BOUNDARY_CONDITIONS_DEPEND_ON_HALO][level];
              prevInd[3] = indexInList[CTaskDefinition::SURFACE_RESIDUAL_HALO_ELEMENTS][level];
              prevInd[4] = indexInList[CTaskDefinition::SURFACE_RESIDUAL_OWNED_ELEMENTS][level];
            }

            /* Create the task. */
            indexInList[CTaskDefinition::ADER_TIME_INTERPOLATE_OWNED_ELEMENTS][level] = (int)tasksList.size();
            tasksList.push_back(CTaskDefinition(CTaskDefinition::ADER_TIME_INTERPOLATE_OWNED_ELEMENTS,
                                                level, prevInd[0], prevInd[1], prevInd[2], prevInd[3], prevInd[4]));

            /* The info on the integration point and whether or not this
               time integration corresponds to the second part for the
               adjacent elements must be added to the task just created. */
            tasksList.back().intPointADER          = intPoint;
            tasksList.back().secondPartTimeIntADER = level < timeLevel;

            /* If artificial viscosity is used for the shock capturing
               terms, these terms must be computed for the owned elements,
               including the adjacent ones. */
            prevInd[0] = indexInList[CTaskDefinition::ADER_TIME_INTERPOLATE_OWNED_ELEMENTS][level];
            indexInList[CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_OWNED_ELEMENTS][level] = (int)tasksList.size();
            tasksList.push_back(CTaskDefinition(CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_OWNED_ELEMENTS,
                                                level, prevInd[0]));
          }

          /* The solution of the halo elements of the current time level and
             the adjacent halo elements of the next time level must be
             interpolated to the integration point intPoint. Check if these
             elements are present. */
          if( interpolHaloElem[level] ) {

            unsigned long nAdjHaloElem = 0;
            if(level < (nTimeLevels-1))
              nAdjHaloElem = haloElemAdjLowTimeLevel[level+1].size();

            unsigned long nHaloElem = nVolElemHaloPerTimeLevel[level+1]
                                    - nVolElemHaloPerTimeLevel[level];

            /* Determine the tasks that should have been completed before this
               task can be carried out. This depends on the time integration
               point. */
            if(intPoint == 0) {

              /* First time integration point. The communication of the
                 halo elements must have been completed before this task
                 can be carried out. */
             if( nHaloElem )
               prevInd[0] = indexInList[CTaskDefinition::COMPLETE_MPI_COMMUNICATION][level];
             else
               prevInd[0] = -1;

             if( nAdjHaloElem )
               prevInd[1] = indexInList[CTaskDefinition::COMPLETE_MPI_COMMUNICATION][level+1];
             else
               prevInd[1] = -1;
            }
            else {

              /* Not the first integration point. The residual computation of
                 the previous integration point must have been completed, because
                 the solution in the work vectors is overwritten. */
              prevInd[0] = indexInList[CTaskDefinition::SURFACE_RESIDUAL_HALO_ELEMENTS][level];
              prevInd[1] = indexInList[CTaskDefinition::BOUNDARY_CONDITIONS_DEPEND_ON_HALO][level];
            }

            /* Create the task. */
            indexInList[CTaskDefinition::ADER_TIME_INTERPOLATE_HALO_ELEMENTS][level] = (int)tasksList.size();
            tasksList.push_back(CTaskDefinition(CTaskDefinition::ADER_TIME_INTERPOLATE_HALO_ELEMENTS,
                                                level, prevInd[0], prevInd[1]));

            /* The info on the integration point and whether or not this
               time integration corresponds to the second part for the
               adjacent elements must be added to the task just created. */
            tasksList.back().intPointADER          = intPoint;
            tasksList.back().secondPartTimeIntADER = level < timeLevel;

            /* If artificial viscosity is used for the shock capturing
               terms, these terms must be computed for the halo elements,
               including the adjacent ones. */
            prevInd[0] = indexInList[CTaskDefinition::ADER_TIME_INTERPOLATE_HALO_ELEMENTS][level];
            indexInList[CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_HALO_ELEMENTS][level] = (int)tasksList.size();
            tasksList.push_back(CTaskDefinition(CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_HALO_ELEMENTS,
                                                level, prevInd[0]));
          }

          /* Check whether there are boundary conditions that involve halo elements. */
          if( BCDependOnHalos[level] ) {

            /* Create the dependency list for this task. */
            prevInd[0] = indexInList[CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_OWNED_ELEMENTS][level];
            prevInd[1] = indexInList[CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_HALO_ELEMENTS][level];

            /* Create the task for the boundary conditions that involve halo elements. */
            indexInList[CTaskDefinition::BOUNDARY_CONDITIONS_DEPEND_ON_HALO][level] = (int)tasksList.size();
            tasksList.push_back(CTaskDefinition(CTaskDefinition::BOUNDARY_CONDITIONS_DEPEND_ON_HALO,
                                                level, prevInd[0], prevInd[1]));
          }

          /* Compute the surface residuals for this time level that involve
             halo elements, if present. Afterwards, accumulate the space time
             residuals for the halo elements. */
          if(nMatchingInternalFacesWithHaloElem[level+1] >
             nMatchingInternalFacesWithHaloElem[level]) {

            /* Create the dependencies for the surface residual part that involve
               halo elements. For all but the first integration point, make sure
               that the previous residual is already accumulated, because this
               task will overwrite that residual. */
            prevInd[0] = indexInList[CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_OWNED_ELEMENTS][level];
            prevInd[1] = indexInList[CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_HALO_ELEMENTS][level];

            if(intPoint == 0)
              prevInd[2] = -1;
            else
              prevInd[2] = indexInList[CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_HALO_ELEMENTS][level];

            /* Create the task for the surface residual. */
            indexInList[CTaskDefinition::SURFACE_RESIDUAL_HALO_ELEMENTS][level] = (int)tasksList.size();
            tasksList.push_back(CTaskDefinition(CTaskDefinition::SURFACE_RESIDUAL_HALO_ELEMENTS,
                                                level, prevInd[0], prevInd[1], prevInd[2]));

            /* Create the task to accumulate the surface residuals of the halo
               elements. Make sure to set the integration point for this task. */
            prevInd[0] = indexInList[CTaskDefinition::SURFACE_RESIDUAL_HALO_ELEMENTS][level];
            indexInList[CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_HALO_ELEMENTS][level] = (int)tasksList.size();
            tasksList.push_back(CTaskDefinition(CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_HALO_ELEMENTS,
                                                level, prevInd[0]));
            tasksList.back().intPointADER = intPoint;
          }

          /* If this is the last integration point, initiate the reverse
             communication for the residuals of this time level, if there
             is something to communicate. */
          if(intPoint == (nTimeIntegrationPoints-1)) {

            /* Check if there is something to communicate. Despite the fact
               that self communication takes place when the reverse communication
               is completed, a check for self communication is still needed,
               because of the periodic transformations. */
            bool commData = false;

#ifdef HAVE_MPI
            if( commRequests[level].size() ) commData = true;
#endif
            if(commData || elementsSendSelfComm[level].size()) {

              /* Determine the dependencies before the reverse communication
                 can start. */
              prevInd[0] = indexInList[CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_HALO_ELEMENTS][level];
              if( haloElemAdjLowTimeLevel[level].size() )
                prevInd[1] = indexInList[CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_HALO_ELEMENTS][level-1];
              else
                prevInd[1] = -1;

              prevInd[2] = indexInList[CTaskDefinition::COMPLETE_MPI_COMMUNICATION][level];

              /* Create the task. */
              indexInList[CTaskDefinition::INITIATE_REVERSE_MPI_COMMUNICATION][level] = (int)tasksList.size();
              tasksList.push_back(CTaskDefinition(CTaskDefinition::INITIATE_REVERSE_MPI_COMMUNICATION,
                                                  level, prevInd[0], prevInd[1], prevInd[2]));
            }
          }

          /* Compute the contribution of the volume integral, boundary conditions that only
             involve owned elements and surface integral between owned elements for this
             time level, if present. */
          if( nOwnedElem ) {

            /* Create the dependencies for this task. The shock capturing viscosity of
               the owned elements must be completed and for all integration points but
               the first, the computed residuals must have been accumulated in the space
               time residual, because the tasks below will overwrite the residuals. */
            prevInd[0] = indexInList[CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_OWNED_ELEMENTS][level];
            if(intPoint == 0)
              prevInd[1] = -1;
            else
              prevInd[1] = indexInList[CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_OWNED_ELEMENTS][level];

            /* Create the tasks. */
            indexInList[CTaskDefinition::VOLUME_RESIDUAL][level] = (int)tasksList.size();
            tasksList.push_back(CTaskDefinition(CTaskDefinition::VOLUME_RESIDUAL, level,
                                                prevInd[0], prevInd[1]));

            indexInList[CTaskDefinition::BOUNDARY_CONDITIONS_DEPEND_ON_OWNED][level] = (int)tasksList.size();
            tasksList.push_back(CTaskDefinition(CTaskDefinition::BOUNDARY_CONDITIONS_DEPEND_ON_OWNED, level,
                                                prevInd[0], prevInd[1]));

            if(nMatchingInternalFacesLocalElem[level+1] >
               nMatchingInternalFacesLocalElem[level]) {
              indexInList[CTaskDefinition::SURFACE_RESIDUAL_OWNED_ELEMENTS][level] = (int)tasksList.size();
              tasksList.push_back(CTaskDefinition(CTaskDefinition::SURFACE_RESIDUAL_OWNED_ELEMENTS,
                                                  level, prevInd[0], prevInd[1]));
            }
          }

          /* Accumulate the space time residuals for the owned elements of this
             level, if these elements are present. */
          if(nAdjOwnedElem || nOwnedElem) {

            /* Create the dependencies for this task. */
            prevInd[0] = indexInList[CTaskDefinition::VOLUME_RESIDUAL][level];
            prevInd[1] = indexInList[CTaskDefinition::BOUNDARY_CONDITIONS_DEPEND_ON_OWNED][level];
            prevInd[1] = indexInList[CTaskDefinition::BOUNDARY_CONDITIONS_DEPEND_ON_HALO][level];
            prevInd[3] = indexInList[CTaskDefinition::SURFACE_RESIDUAL_OWNED_ELEMENTS][level];
            prevInd[4] = indexInList[CTaskDefinition::SURFACE_RESIDUAL_HALO_ELEMENTS][level];

            /* Create the task. */
            indexInList[CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_OWNED_ELEMENTS][level] = (int)tasksList.size();
            tasksList.push_back(CTaskDefinition(CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_OWNED_ELEMENTS,
                                                level, prevInd[0], prevInd[1], prevInd[2], prevInd[3], prevInd[4]));
            tasksList.back().intPointADER = intPoint;
          }

          /* If this is the last integration point, complete the reverse
             communication for the residuals of this time level, if there
             is something to communicate. */
          if(intPoint == (nTimeIntegrationPoints-1)) {

            bool commData = false;

#ifdef HAVE_MPI
            if( commRequests[level].size() ) commData = true;
#endif
            if(commData || elementsSendSelfComm[level].size()) {

              /* Determine the dependencies before the reverse communication
                 can be completed. */
              prevInd[0] = indexInList[CTaskDefinition::INITIATE_REVERSE_MPI_COMMUNICATION][level];
              prevInd[1] = indexInList[CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_OWNED_ELEMENTS][level];
              if( ownedElemAdjLowTimeLevel[level].size() )
                prevInd[2] = indexInList[CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_OWNED_ELEMENTS][level-1];
              else
                prevInd[2] = -1;

              /* Create the task. */
              indexInList[CTaskDefinition::COMPLETE_REVERSE_MPI_COMMUNICATION][level] = (int)tasksList.size();
              tasksList.push_back(CTaskDefinition(CTaskDefinition::COMPLETE_REVERSE_MPI_COMMUNICATION,
                                                  level, prevInd[0], prevInd[1], prevInd[2]));
            }
          }

        } /* End loop time level. */

      } /* End loop time integration points. */

      /* Loop again over the number of active time levels for this subtime
         step and compute the update for the DOFs of the owned elements. */
      for(unsigned short level=0; level<=timeLevel; ++level) {
        if(nVolElemOwnedPerTimeLevel[level+1] > nVolElemOwnedPerTimeLevel[level]) {

          /* Updates must be carried out for this time level. First multiply the
             residuals by the inverse of the mass matrix. This task can be
             completed if the communication of this level has been completed.
             However, it is possible that no communication is needed and
             therefore the accumulation of the space time residuals is also
             added to the dependency list. */
          prevInd[0] = indexInList[CTaskDefinition::COMPLETE_REVERSE_MPI_COMMUNICATION][level];
          prevInd[1] = indexInList[CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_OWNED_ELEMENTS][level];
          if( ownedElemAdjLowTimeLevel[level].size() )
            prevInd[2] = indexInList[CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_OWNED_ELEMENTS][level-1];
          else
            prevInd[2] = -1;

          indexInList[CTaskDefinition::MULTIPLY_INVERSE_MASS_MATRIX][level] = (int)tasksList.size();
          tasksList.push_back(CTaskDefinition(CTaskDefinition::MULTIPLY_INVERSE_MASS_MATRIX,
                                               level, prevInd[0], prevInd[1], prevInd[2]));

          /* Compute the new state vector for this time level. */
          prevInd[0] = indexInList[CTaskDefinition::MULTIPLY_INVERSE_MASS_MATRIX][level];
          indexInList[CTaskDefinition::ADER_UPDATE_SOLUTION][level] = (int)tasksList.size();
          tasksList.push_back(CTaskDefinition(CTaskDefinition::ADER_UPDATE_SOLUTION,
                                               level, prevInd[0]));
        }
      }

    } /* End loop subtime steps. */


    // EXTRA FOR DEBUGGING
/*  for(int i=0; i<size; ++i) {
      if(i == rank) {
        cout << endl;
        cout << "Task list for rank " << rank << endl;
        cout << "------------------------------------------------" << endl;
        cout << "Number of tasks: " << tasksList.size() << endl;
        for(unsigned long j=0; j<tasksList.size(); ++j) {

          cout << "Task " << j << ": ";
          switch( tasksList[j].task ) {
            case CTaskDefinition::NO_TASK: cout << "NO_TASK" << endl; break;
            case CTaskDefinition::ADER_PREDICTOR_STEP_COMM_ELEMENTS: cout << "ADER_PREDICTOR_STEP_COMM_ELEMENTS" << endl; break;
            case CTaskDefinition::ADER_PREDICTOR_STEP_INTERNAL_ELEMENTS: cout << "ADER_PREDICTOR_STEP_INTERNAL_ELEMENTS" << endl; break;
            case CTaskDefinition::INITIATE_MPI_COMMUNICATION: cout << "INITIATE_MPI_COMMUNICATION" << endl; break;
            case CTaskDefinition::COMPLETE_MPI_COMMUNICATION: cout << "COMPLETE_MPI_COMMUNICATION" << endl; break;
            case CTaskDefinition::INITIATE_REVERSE_MPI_COMMUNICATION: cout << "INITIATE_REVERSE_MPI_COMMUNICATION" << endl; break;
            case CTaskDefinition::COMPLETE_REVERSE_MPI_COMMUNICATION: cout << "COMPLETE_REVERSE_MPI_COMMUNICATION" << endl; break;
            case CTaskDefinition::ADER_TIME_INTERPOLATE_OWNED_ELEMENTS: cout << "ADER_TIME_INTERPOLATE_OWNED_ELEMENTS" << endl; break;
            case CTaskDefinition::ADER_TIME_INTERPOLATE_HALO_ELEMENTS: cout << "ADER_TIME_INTERPOLATE_HALO_ELEMENTS" << endl; break;
            case CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_OWNED_ELEMENTS: cout << "SHOCK_CAPTURING_VISCOSITY_OWNED_ELEMENTS" << endl; break;
            case CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_HALO_ELEMENTS: cout << "SHOCK_CAPTURING_VISCOSITY_HALO_ELEMENTS" << endl; break;
            case CTaskDefinition::VOLUME_RESIDUAL: cout << "VOLUME_RESIDUAL" << endl; break;
            case CTaskDefinition::SURFACE_RESIDUAL_OWNED_ELEMENTS: cout << "SURFACE_RESIDUAL_OWNED_ELEMENTS" << endl; break;
            case CTaskDefinition::SURFACE_RESIDUAL_HALO_ELEMENTS: cout << "SURFACE_RESIDUAL_HALO_ELEMENTS" << endl; break;
            case CTaskDefinition::BOUNDARY_CONDITIONS_DEPEND_ON_OWNED: cout << "BOUNDARY_CONDITIONS_DEPEND_ON_OWNED" << endl; break;
            case CTaskDefinition::BOUNDARY_CONDITIONS_DEPEND_ON_HALO: cout << "BOUNDARY_CONDITIONS_DEPEND_ON_HALO" << endl; break;
            case CTaskDefinition::SUM_UP_RESIDUAL_CONTRIBUTIONS_OWNED_ELEMENTS: cout << "SUM_UP_RESIDUAL_CONTRIBUTIONS_OWNED_ELEMENTS" << endl; break;
            case CTaskDefinition::SUM_UP_RESIDUAL_CONTRIBUTIONS_HALO_ELEMENTS: cout << "SUM_UP_RESIDUAL_CONTRIBUTIONS_HALO_ELEMENTS" << endl; break;
            case CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_OWNED_ELEMENTS: cout << "ADER_ACCUMULATE_SPACETIME_RESIDUAL_OWNED_ELEMENTS" << endl; break;
            case CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_HALO_ELEMENTS: cout << "ADER_ACCUMULATE_SPACETIME_RESIDUAL_HALO_ELEMENTS" << endl; break;
            case CTaskDefinition::MULTIPLY_INVERSE_MASS_MATRIX: cout << "MULTIPLY_INVERSE_MASS_MATRIX" << endl; break;
            case CTaskDefinition::ADER_UPDATE_SOLUTION: cout << "ADER_UPDATE_SOLUTION" << endl; break;
            default: cout << "This cannot happen" << endl;
          }
          cout << " Time level: " << tasksList[j].timeLevel
               << " Integration point: " << tasksList[j].intPointADER
               << " Second part: " << tasksList[j].secondPartTimeIntADER << endl;
          cout << " Depends on tasks:";
          for(unsigned short k=0; k<tasksList[j].nIndMustBeCompleted; ++k)
            cout << " " << tasksList[j].indMustBeCompleted[k];
          cout << endl << endl;
        }

      }

#ifdef HAVE_MPI
      SU2_MPI::Barrier(SU2_MPI::GetComm());
#endif
    }

    cout << "CFEM_DG_EulerSolver::SetUpTaskList: ADER tasklist printed" << endl;
    exit(1); */

    // END EXTRA FOR DEBUGGING.
  }
  else {

    /*------------------------------------------------------------------------*/
    /* Standard time integration scheme for which the spatial residual must   */
    /* be computed for the DOFS of the owned elements. This results in a      */
    /* relatively short tasks list, which can be set easily.                  */
    /*------------------------------------------------------------------------*/

    tasksList.push_back(CTaskDefinition(CTaskDefinition::INITIATE_MPI_COMMUNICATION,                   0));                   // Task  0
    tasksList.push_back(CTaskDefinition(CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_OWNED_ELEMENTS,     0));                   // Task  1
    tasksList.push_back(CTaskDefinition(CTaskDefinition::VOLUME_RESIDUAL,                              0,  1));               // Task  2
    tasksList.push_back(CTaskDefinition(CTaskDefinition::COMPLETE_MPI_COMMUNICATION,                   0,  0));               // Task  3
    tasksList.push_back(CTaskDefinition(CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_HALO_ELEMENTS,      0,  3));               // Task  4
    tasksList.push_back(CTaskDefinition(CTaskDefinition::SURFACE_RESIDUAL_HALO_ELEMENTS,               0,  1, 4));            // Task  5
    tasksList.push_back(CTaskDefinition(CTaskDefinition::BOUNDARY_CONDITIONS_DEPEND_ON_HALO,           0,  1, 3));            // Task  6
    tasksList.push_back(CTaskDefinition(CTaskDefinition::SUM_UP_RESIDUAL_CONTRIBUTIONS_HALO_ELEMENTS,  0,  5));               // Task  7
    tasksList.push_back(CTaskDefinition(CTaskDefinition::INITIATE_REVERSE_MPI_COMMUNICATION,           0,  7));               // Task  8
    tasksList.push_back(CTaskDefinition(CTaskDefinition::SURFACE_RESIDUAL_OWNED_ELEMENTS,              0,  1));               // Task  9
    tasksList.push_back(CTaskDefinition(CTaskDefinition::BOUNDARY_CONDITIONS_DEPEND_ON_OWNED,          0,  1));               // Task 10
    tasksList.push_back(CTaskDefinition(CTaskDefinition::COMPLETE_REVERSE_MPI_COMMUNICATION,           0,  2, 8));            // Task 11
    tasksList.push_back(CTaskDefinition(CTaskDefinition::SUM_UP_RESIDUAL_CONTRIBUTIONS_OWNED_ELEMENTS, 0,  2, 6, 9, 10, 11)); // Task 12
    tasksList.push_back(CTaskDefinition(CTaskDefinition::MULTIPLY_INVERSE_MASS_MATRIX,                 0, 12));               // Task 13
  }
}

void CFEM_DG_EulerSolver::ConservativeToEntropyVariables(ColMajorMatrix<su2double> &sol) {

  /*--- Easier storage of the number of items to be converted and
        the inverse of gamma-1. ---*/
  const unsigned short nItems = sol.rows();
  const su2double ovgm1       = 1.0/Gamma_Minus_One;

  /*--- Make a distinction between two and three space dimensions
        in order to have the most efficient code. ---*/
  switch( nDim ) {

    case 2: {

      /*--- Two dimensional simulation. Loop over the number of entities
            and carry out the conversion. ---*/
      SU2_OMP_SIMD_IF_NOT_AD
      for(unsigned short i=0; i<nItems; ++i) {
        const su2double rho    = sol(i,0);
        const su2double rhoInv = 1.0/rho;
        const su2double u      = rhoInv*sol(i,1);
        const su2double v      = rhoInv*sol(i,2);
        const su2double kin    = 0.5*(u*u + v*v);
        const su2double p      = Gamma_Minus_One*(sol(i,3) - rho*kin);
        const su2double pInv   = 1.0/p;
        const su2double s      = log(p*pow(rhoInv,Gamma));
        const su2double abv    = pInv*rho;

        sol(i,0) =  (Gamma-s)*ovgm1 - abv*kin;
        sol(i,1) =  abv*u;
        sol(i,2) =  abv*v;
        sol(i,3) = -abv;
      }

      break;
    }

    case 3: {

      /*--- Two dimensional simulation. Loop over the number of entities
            and carry out the conversion. ---*/
      SU2_OMP_SIMD_IF_NOT_AD
      for(unsigned short i=0; i<nItems; ++i) {
        const su2double rho    = sol(i,0);
        const su2double rhoInv = 1.0/rho;
        const su2double u      = rhoInv*sol(i,1);
        const su2double v      = rhoInv*sol(i,2);
        const su2double w      = rhoInv*sol(i,3);
        const su2double kin    = 0.5*(u*u + v*v + w*w);
        const su2double p      = Gamma_Minus_One*(sol(i,4) - rho*kin);
        const su2double pInv   = 1.0/p;
        const su2double s      = log(p*pow(rhoInv,Gamma));
        const su2double abv    = pInv*rho;

        sol(i,0) =  (Gamma-s)*ovgm1 - pInv*rho*kin;
        sol(i,1) =  abv*u;
        sol(i,2) =  abv*v;
        sol(i,3) =  abv*w;
        sol(i,4) = -abv;
      }

      break;
    }
  }
}

void CFEM_DG_EulerSolver::Initiate_MPI_Communication(CConfig *config,
                                                     const unsigned short timeLevel) {
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

bool CFEM_DG_EulerSolver::Complete_MPI_Communication(CConfig *config,
                                                     const unsigned short timeLevel,
                                                     const bool commMustBeCompleted) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
  return false;
}

void CFEM_DG_EulerSolver::Initiate_MPI_ReverseCommunication(CConfig *config,
                                                            const unsigned short timeLevel) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

bool CFEM_DG_EulerSolver::Complete_MPI_ReverseCommunication(CConfig *config,
                                                            const unsigned short timeLevel,
                                                            const bool commMustBeCompleted) {
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
  return false;
}

void CFEM_DG_EulerSolver::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long TimeIter) {

  /*--- Check if a verification solution is to be computed. ---*/
  if ((VerificationSolution)  && (TimeIter == 0)) {

    /*--- Start OpenMP parallel region. ---*/
    SU2_OMP_PARALLEL {

      /*--- Loop over the owned elements. ---*/
#ifdef HAVE_OMP
      const size_t omp_chunk_size_vol = computeStaticChunkSize(nVolElemOwned, omp_get_num_threads(), 64);
#endif
      SU2_OMP_FOR_STAT(omp_chunk_size_vol)
      for(unsigned long i=0; i<nVolElemOwned; ++i) {

        /*--- Define the help variables to store the coordinates
              and the solution in the nodal DOFs. ---*/
        su2double coor[3] = {0.0}, solDOF[5] = {0.0};

        /*--- Loop over the DOFs of this element. ---*/
        const unsigned short nDOFs = volElem[i].standardElemFlow->GetNDOFs();
        for(unsigned short j=0; j<nDOFs; ++j) {

          /*--- Store the coordinates of the DOFs in coor. ---*/
          for(unsigned short iDim=0; iDim<nDim; ++iDim)
            coor[iDim] = volElem[i].coorSolDOFs(j,iDim);

          /*--- Compute the initial condition provided by the verification solution
                class. This can be the exact solution, but this is not necessary.
                Note that the conservative variables are computed. ---*/
          VerificationSolution->GetInitialCondition(coor, solDOF);

          /*--- Store the conservative variables in the DOFs for the element
                for the time being. ---*/
          for(unsigned short iVar=0; iVar<nVar; ++iVar)
            volElem[i].solDOFs(j,iVar) = solDOF[iVar];
        }

        /*--- Convert the conservative variables to entropy variables. ---*/
        ConservativeToEntropyVariables(volElem[i].solDOFs);

        /*--- Convert the nodal solution to the modal solution. ---*/
        volElem[i].NodalToModalFlow();
      }
    } // end SU2_OMP_PARALLEL
  }
}

void CFEM_DG_EulerSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iStep, unsigned short RunTime_EqSystem, bool Output) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                      unsigned short iMesh) { }

void CFEM_DG_EulerSolver::ComputeSpatialJacobian(CGeometry *geometry,  CSolver **solver_container,
                                                 CNumerics **numerics, CConfig *config,
                                                 unsigned short iMesh, unsigned short RunTime_EqSystem) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::Set_OldSolution() {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::Set_NewSolution() {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                    unsigned short iMesh, unsigned long Iteration) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::CheckTimeSynchronization(CConfig         *config,
                                                   const su2double TimeSync,
                                                   su2double       &timeEvolved,
                                                   bool            &syncTimeReached) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::ProcessTaskList_DG(CGeometry *geometry,  CSolver **solver_container,
                                             CNumerics **numerics, CConfig *config,
                                             unsigned short iMesh) {
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::ADER_SpaceTimeIntegration(CGeometry *geometry,  CSolver **solver_container,
                                                    CNumerics **numerics, CConfig *config,
                                                    unsigned short iMesh, unsigned short RunTime_EqSystem) {
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::TolerancesADERPredictorStep(void) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::ADER_DG_PredictorStep(CConfig             *config,
                                                const unsigned long elemBeg,
                                                const unsigned long elemEnd,
                                                su2double           *workArray) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::ADER_DG_AliasedPredictorResidual_2D(CConfig              *config,
                                                              CVolumeElementFEM_DG *elem,
                                                              const su2double      *sol,
                                                              const unsigned short nSimul,
                                                              const unsigned short NPad,
                                                              su2double            *res,
                                                              su2double            *work) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::ADER_DG_AliasedPredictorResidual_3D(CConfig              *config,
                                                              CVolumeElementFEM_DG *elem,
                                                              const su2double      *sol,
                                                              const unsigned short nSimul,
                                                              const unsigned short NPad,
                                                              su2double            *res,
                                                              su2double            *work) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::ADER_DG_NonAliasedPredictorResidual_2D(CConfig              *config,
                                                                 CVolumeElementFEM_DG *elem,
                                                                 const su2double      *sol,
                                                                 const unsigned short nSimul,
                                                                 const unsigned short NPad,
                                                                 su2double            *res,
                                                                 su2double            *work) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::ADER_DG_NonAliasedPredictorResidual_3D(CConfig              *config,
                                                                 CVolumeElementFEM_DG *elem,
                                                                 const su2double      *sol,
                                                                 const unsigned short nSimul,
                                                                 const unsigned short NPad,
                                                                 su2double            *res,
                                                                 su2double            *work) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::ADER_DG_TimeInterpolatePredictorSol(CConfig             *config,
                                                              const unsigned short iTime,
                                                              const unsigned long  elemBeg,
                                                              const unsigned long  elemEnd,
                                                              const unsigned long  nAdjElem,
                                                              const unsigned long  *adjElem,
                                                              const bool           secondPartTimeInt,
                                                              su2double            *solTimeLevel) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::Shock_Capturing_DG(CConfig             *config,
                                             const unsigned long elemBeg,
                                             const unsigned long elemEnd,
                                             su2double           *workArray) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::Volume_Residual(CConfig             *config,
                                          const unsigned long elemBeg,
                                          const unsigned long elemEnd,
                                          su2double           *workArray) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::Boundary_Conditions(const unsigned short timeLevel,
                                              CConfig              *config,
                                              CNumerics            **numerics,
                                              const bool           haloInfoNeededForBC,
                                              su2double            *workArray){

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::ResidualFaces(CConfig             *config,
                                        const unsigned long indFaceBeg,
                                        const unsigned long indFaceEnd,
                                        unsigned long       &indResFaces,
                                        CNumerics           *numerics,
                                        su2double           *workArray) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::InviscidFluxesInternalMatchingFace(
                                              CConfig              *config,
                                              const unsigned long  lBeg,
                                              const unsigned long  lEnd,
                                              const unsigned short NPad,
                                              su2double            *solIntL,
                                              su2double            *solIntR,
                                              su2double            *fluxes,
                                              CNumerics            *numerics) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::AccumulateSpaceTimeResidualADEROwnedElem(
                                                     CConfig             *config,
                                                     const unsigned short timeLevel,
                                                     const unsigned short intPoint) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::AccumulateSpaceTimeResidualADERHaloElem(
                                                     CConfig             *config,
                                                     const unsigned short timeLevel,
                                                     const unsigned short intPoint) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::CreateFinalResidual(const unsigned short timeLevel,
                                              const bool ownedElements) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::MultiplyResidualByInverseMassMatrix(
                                              CConfig            *config,
                                              const bool          useADER,
                                              const unsigned long elemBeg,
                                              const unsigned long elemEnd,
                                              su2double           *workArray) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::Pressure_Forces(const CGeometry *geometry, const CConfig *config) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::ExplicitRK_Iteration(CGeometry *geometry, CSolver **solver_container,
                                               CConfig *config, unsigned short iRKStep) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::ClassicalRK4_Iteration(CGeometry *geometry, CSolver **solver_container,
                                                 CConfig *config, unsigned short iRKStep) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::SetResidual_RMS_FEM(CGeometry *geometry,
                                              CConfig *config) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::ComputeVerificationError(CGeometry *geometry,
                                                   CConfig   *config) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::ADER_DG_Iteration(const unsigned long elemBeg,
                                            const unsigned long elemEnd) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::BoundaryStates_Euler_Wall(CConfig                  *config,
                                                    const unsigned short     nFaceSimul,
                                                    const unsigned short     NPad,
                                                    const CSurfaceElementFEM *surfElem,
                                                    const su2double          *solIntL,
                                                    su2double                *solIntR) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::BoundaryStates_Inlet(CConfig                  *config,
                                               const unsigned short     nFaceSimul,
                                               const unsigned short     NPad,
                                               const CSurfaceElementFEM *surfElem,
                                               unsigned short           val_marker,
                                               const su2double          *solIntL,
                                               su2double                *solIntR) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::BoundaryStates_Outlet(CConfig                  *config,
                                                const unsigned short     nFaceSimul,
                                                const unsigned short     NPad,
                                                const CSurfaceElementFEM *surfElem,
                                                unsigned short           val_marker,
                                                const su2double          *solIntL,
                                                su2double                *solIntR) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::BoundaryStates_Riemann(CConfig                  *config,
                                                 const unsigned short     nFaceSimul,
                                                 const unsigned short     NPad,
                                                 const CSurfaceElementFEM *surfElem,
                                                 unsigned short           val_marker,
                                                 const su2double          *solIntL,
                                                 su2double                *solIntR) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::BC_Euler_Wall(CConfig                  *config,
                                        const unsigned long      surfElemBeg,
                                        const unsigned long      surfElemEnd,
                                        const CSurfaceElementFEM *surfElem,
                                        su2double                *resFaces,
                                        CNumerics                *conv_numerics,
                                        su2double                *workArray) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::BC_Far_Field(CConfig                  *config,
                                       const unsigned long      surfElemBeg,
                                       const unsigned long      surfElemEnd,
                                       const CSurfaceElementFEM *surfElem,
                                       su2double                *resFaces,
                                       CNumerics                *conv_numerics,
                                       su2double                *workArray) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::BC_Sym_Plane(CConfig                  *config,
                                       const unsigned long      surfElemBeg,
                                       const unsigned long      surfElemEnd,
                                       const CSurfaceElementFEM *surfElem,
                                       su2double                *resFaces,
                                       CNumerics                *conv_numerics,
                                       su2double                *workArray) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::BC_Supersonic_Outlet(CConfig                  *config,
                                               const unsigned long      surfElemBeg,
                                               const unsigned long      surfElemEnd,
                                               const CSurfaceElementFEM *surfElem,
                                               su2double                *resFaces,
                                               CNumerics                *conv_numerics,
                                               su2double                *workArray) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::BC_Inlet(CConfig                  *config,
                                   const unsigned long      surfElemBeg,
                                   const unsigned long      surfElemEnd,
                                   const CSurfaceElementFEM *surfElem,
                                   su2double                *resFaces,
                                   CNumerics                *conv_numerics,
                                   unsigned short           val_marker,
                                   su2double                *workArray) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::BC_Outlet(CConfig                  *config,
                                    const unsigned long      surfElemBeg,
                                    const unsigned long      surfElemEnd,
                                    const CSurfaceElementFEM *surfElem,
                                    su2double                *resFaces,
                                    CNumerics                *conv_numerics,
                                    unsigned short           val_marker,
                                    su2double                *workArray) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::BC_Riemann(CConfig                  *config,
                                     const unsigned long      surfElemBeg,
                                     const unsigned long      surfElemEnd,
                                     const CSurfaceElementFEM *surfElem,
                                     su2double                *resFaces,
                                     CNumerics                *conv_numerics,
                                     unsigned short           val_marker,
                                     su2double                *workArray) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::BC_Custom(CConfig                  *config,
                                    const unsigned long      surfElemBeg,
                                    const unsigned long      surfElemEnd,
                                    const CSurfaceElementFEM *surfElem,
                                    su2double                *resFaces,
                                    CNumerics                *conv_numerics,
                                    su2double                *workArray) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::ResidualInviscidBoundaryFace(
                                      CConfig                  *config,
                                      const unsigned short     nFaceSimul,
                                      const unsigned short     NPad,
                                      CNumerics                *conv_numerics,
                                      const CSurfaceElementFEM *surfElem,
                                      su2double                *solInt0,
                                      su2double                *solInt1,
                                      su2double                *fluxes,
                                      su2double                *resFaces,
                                      unsigned long            &indResFaces) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::LeftStatesIntegrationPointsBoundaryFace(
                                             CConfig                  *config,
                                             const unsigned short     nFaceSimul,
                                             const unsigned short     NPad,
                                             const CSurfaceElementFEM *surfElem,
                                             su2double                *solFace,
                                             su2double                *solIntL) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::ComputeInviscidFluxesFace(CConfig              *config,
                                                    const unsigned short nFaceSimul,
                                                    const unsigned short NPad,
                                                    const unsigned long  nPoints,
                                                    const su2double      *normalsFace[],
                                                    const su2double      *gridVelsFace[],
                                                    const su2double      *solL,
                                                    const su2double      *solR,
                                                    su2double            *fluxes,
                                                    CNumerics            *numerics) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::LoadRestart(CGeometry **geometry,
                                      CSolver   ***solver,
                                      CConfig   *config,
                                      int       val_iter,
                                      bool      val_update_geo) {

  /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/
  string restart_filename = config->GetSolution_FileName();
  restart_filename = config->GetFilename(restart_filename, "", val_iter);

  if (config->GetRead_Binary_Restart()) {
    Read_SU2_Restart_Binary(geometry[MESH_0], config, restart_filename);
  } else {
    Read_SU2_Restart_ASCII(geometry[MESH_0], config, restart_filename);
  }

  /*--- Determine a map from the ID of the global DOF to the local
        element and DOF of the element. ---*/
  map<unsigned long, CUnsignedLong2T>  mapGlobalDOFToLocal;
  for(unsigned long i=0; i<nVolElemOwned; ++i) {
    for(unsigned short j=0; j<volElem[i].nDOFsSol; ++j) {
      unsigned long GlobalDOF = volElem[i].offsetDOFsSolGlobal+j;
      mapGlobalDOFToLocal[GlobalDOF] = CUnsignedLong2T(i,j);
    }
  }

  /*--- Skip coordinates ---*/
  unsigned short skipVars = nDim;

  /*--- Loop over the global DOFs and determine whether they are
        stored locally. ---*/
  unsigned long counter = 0;
  for (unsigned long iPoint_Global=0; iPoint_Global<geometry[MESH_0]->GetGlobal_nPointDomain(); ++iPoint_Global) {

    /*--- Check if the DOF is stored on this rank. ---*/
    auto it = mapGlobalDOFToLocal.find(iPoint_Global);
    if(it != mapGlobalDOFToLocal.cend()) {

      /*--- Retrieve the element and DOF inside the element. ---*/
      const unsigned long i = it->second.long0;
      const unsigned long j = it->second.long1;

      /*--- Determine the start index in the read buffer and
            update the counter afterwards. ---*/
      const unsigned long index = counter*Restart_Vars[1] + skipVars;
      ++counter;

      /*--- Store the conservative variables in the local data structures. ---*/
      for(unsigned short iVar=0; iVar<nVar; ++iVar)
        volElem[i].solDOFs(j,iVar) = Restart_Data[index+iVar];
    }
  }

  /*--- Delete the class memory that is used to load the restart. ---*/
  delete[] Restart_Vars;
  delete[] Restart_Data;
  Restart_Vars = nullptr; Restart_Data = nullptr;

  /*--- Detect a wrong solution file ---*/
  unsigned short rbuf_NotMatching = 0;
  if(counter < nDOFsLocOwned) rbuf_NotMatching = 1;

#ifdef HAVE_MPI
  unsigned short sbuf_NotMatching = rbuf_NotMatching;
  SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_MAX, SU2_MPI::GetComm());
#endif

  if (rbuf_NotMatching != 0)
    SU2_MPI::Error(string("The solution file ") + restart_filename.data() +
                   string(" doesn't match with the mesh file!\n") +
                   string("It could be empty lines at the end of the file."),
                   CURRENT_FUNCTION);

  /*--- Start of the parallel region. ---*/
  unsigned long nBadDOFs = 0;
  SU2_OMP_PARALLEL
  {
    /*--- Loop over the owned elements. ---*/
#ifdef HAVE_OMP
    const size_t omp_chunk_size_vol = computeStaticChunkSize(nVolElemOwned, omp_get_num_threads(), 64);
#endif
    SU2_OMP_FOR_STAT_(omp_chunk_size_vol, reduction(+: nBadDOFs))
    for(unsigned long i=0; i<nVolElemOwned; ++i) {

      /*--- Loop over the DOFs of this element, which currently
            stores the conservative variables. ---*/
      for(unsigned short j=0; j<volElem[i].nDOFsSol; ++j) {

        /*--- Compute the primitive variables. ---*/
        const su2double rho    = volElem[i].solDOFs(j,0);
        const su2double rhoInv = 1.0/rho;

        su2double Velocity2 = 0.0;
        for(unsigned short iDim=1; iDim<=nDim; ++iDim) {
          const su2double vel = volElem[i].solDOFs(j,iDim)*rhoInv;
          Velocity2 += vel*vel;
        }

        const su2double StaticEnergy = volElem[i].solDOFs(j,nDim+1)*rhoInv - 0.5*Velocity2;

        GetFluidModel()->SetTDState_rhoe(rho, StaticEnergy);
        const su2double Pressure = GetFluidModel()->GetPressure();
        const su2double Temperature = GetFluidModel()->GetTemperature();

        /*--- Check for negative pressure, density or temperature. ---*/
        if((Pressure < 0.0) || (rho < 0.0) || (Temperature < 0.0)) {

          /*--- Reset the state to the infinity state and update nBadDOFs. ---*/
          for(unsigned short k=0; k<nVar; ++k)
            volElem[i].solDOFs(j,k) = ConsVarFreeStream[k];
          ++nBadDOFs;
        }
      }

      /*--- Convert the conservative variables to entropy variables. ---*/
      ConservativeToEntropyVariables(volElem[i].solDOFs);

      /*--- Convert the nodal solution to the modal solution. ---*/
      volElem[i].NodalToModalFlow();
    }
  } // end SU2_OMP_PARALLEL

  /*--- Warning message about non-physical points ---*/
  if (config->GetComm_Level() == COMM_FULL) {
#ifdef HAVE_MPI
    unsigned long nBadDOFsLoc = nBadDOFs;
    SU2_MPI::Reduce(&nBadDOFsLoc, &nBadDOFs, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, SU2_MPI::GetComm());
#endif

    if((rank == MASTER_NODE) && (nBadDOFs != 0))
      cout << "Warning. The initial solution contains "<< nBadDOFs << " DOFs that are not physical." << endl;
  }
}
