/***********************************************************************
/
/ Example executable using libgrackle
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

#include <grackle.h>

#define mh     1.67262171e-24   
#define kboltz 1.3806504e-16

void malloc_fields(grackle_field_data *my_fields, int field_size, int dust_flag)
{
  my_fields->grid_dimension  = malloc(field_size * sizeof(gr_float));
  my_fields->grid_start      = malloc(field_size * sizeof(gr_float));
  my_fields->grid_end        = malloc(field_size * sizeof(gr_float));

  my_fields->density         = malloc(field_size * sizeof(gr_float));
  my_fields->internal_energy = malloc(field_size * sizeof(gr_float));
  my_fields->x_velocity      = malloc(field_size * sizeof(gr_float));
  my_fields->y_velocity      = malloc(field_size * sizeof(gr_float));
  my_fields->z_velocity      = malloc(field_size * sizeof(gr_float));
  // for primordial_chemistry >= 1
  my_fields->HI_density      = malloc(field_size * sizeof(gr_float));
  my_fields->HII_density     = malloc(field_size * sizeof(gr_float));
  my_fields->HeI_density     = malloc(field_size * sizeof(gr_float));
  my_fields->HeII_density    = malloc(field_size * sizeof(gr_float));
  my_fields->HeIII_density   = malloc(field_size * sizeof(gr_float));
  my_fields->e_density       = malloc(field_size * sizeof(gr_float));
  // for primordial_chemistry >= 2
  my_fields->HM_density      = malloc(field_size * sizeof(gr_float));
  my_fields->H2I_density     = malloc(field_size * sizeof(gr_float));
  my_fields->H2II_density    = malloc(field_size * sizeof(gr_float));
  // for primordial_chemistry >= 3
  my_fields->DI_density      = malloc(field_size * sizeof(gr_float));
  my_fields->DII_density     = malloc(field_size * sizeof(gr_float));
  my_fields->HDI_density     = malloc(field_size * sizeof(gr_float));
  // for metal_cooling = 1
  my_fields->metal_density   = malloc(field_size * sizeof(gr_float));

  // volumetric heating rate (provide in units [erg s^-1 cm^-3])
  my_fields->volumetric_heating_rate = malloc(field_size * sizeof(gr_float));
  // specific heating rate (provide in units [egs s^-1 g^-1]
  my_fields->specific_heating_rate = malloc(field_size * sizeof(gr_float));
  my_fields->temperature_floor = malloc(field_size * sizeof(gr_float));

  // radiative transfer ionization / dissociation rate fields (provide in units [1/s])
  my_fields->RT_HI_ionization_rate = malloc(field_size * sizeof(gr_float));
  my_fields->RT_HeI_ionization_rate = malloc(field_size * sizeof(gr_float));
  my_fields->RT_HeII_ionization_rate = malloc(field_size * sizeof(gr_float));
  my_fields->RT_H2_dissociation_rate = malloc(field_size * sizeof(gr_float));
  // radiative transfer heating rate field (provide in units [erg s^-1 cm^-3])
  my_fields->RT_heating_rate = malloc(field_size * sizeof(gr_float));

  // H2 model
  my_fields->H2_self_shielding_length = malloc(field_size * sizeof(gr_float));
  my_fields->H2_custom_shielding_factor = malloc(field_size * sizeof(gr_float));
  my_fields->isrf_habing = malloc(field_size * sizeof(gr_float));

  if (dust_flag) {
      my_fields->dust_density = malloc(field_size * sizeof(gr_float));
      my_fields->He_gas_metalDensity = malloc(field_size * sizeof(gr_float));
      my_fields->C_gas_metalDensity = malloc(field_size * sizeof(gr_float));
      my_fields->N_gas_metalDensity = malloc(field_size * sizeof(gr_float));
      my_fields->O_gas_metalDensity = malloc(field_size * sizeof(gr_float));
      my_fields->Ne_gas_metalDensity = malloc(field_size * sizeof(gr_float));
      my_fields->Mg_gas_metalDensity = malloc(field_size * sizeof(gr_float));
      my_fields->Si_gas_metalDensity = malloc(field_size * sizeof(gr_float));
      my_fields->S_gas_metalDensity = malloc(field_size * sizeof(gr_float));
      my_fields->Ca_gas_metalDensity = malloc(field_size * sizeof(gr_float));
      my_fields->Fe_gas_metalDensity = malloc(field_size * sizeof(gr_float));
      my_fields->He_dust_metalDensity = malloc(field_size * sizeof(gr_float));
      my_fields->C_dust_metalDensity = malloc(field_size * sizeof(gr_float));
      my_fields->N_dust_metalDensity = malloc(field_size * sizeof(gr_float));
      my_fields->O_dust_metalDensity = malloc(field_size * sizeof(gr_float));
      my_fields->Ne_dust_metalDensity = malloc(field_size * sizeof(gr_float));
      my_fields->Mg_dust_metalDensity = malloc(field_size * sizeof(gr_float));
      my_fields->Si_dust_metalDensity = malloc(field_size * sizeof(gr_float));
      my_fields->S_dust_metalDensity = malloc(field_size * sizeof(gr_float));
      my_fields->Ca_dust_metalDensity = malloc(field_size * sizeof(gr_float));
      my_fields->Fe_dust_metalDensity = malloc(field_size * sizeof(gr_float));
      my_fields->SNe_ThisTimeStep = malloc(field_size * sizeof(gr_float));
  }
  return;
}

void init_fields(grackle_field_data *my_fields, chemistry_data *grackle_data, code_units my_units, int field_size, double temperature_units, double initial_redshift, double dt) 
{
  int i;
  double tiny_number = 1.e-20;

  my_fields->grid_rank = 3;
  my_fields->grid_dimension = malloc(my_fields->grid_rank * sizeof(int));
  my_fields->grid_start = malloc(my_fields->grid_rank * sizeof(int));
  my_fields->grid_end = malloc(my_fields->grid_rank * sizeof(int));
  my_fields->grid_dx = 0.0; // used only for H2 self-shielding approximation

  for (i = 0;i < 3;i++) {
    my_fields->grid_dimension[i] = 1; // the active dimension not including ghost zones.
    my_fields->grid_start[i] = 0;
    my_fields->grid_end[i] = 0;
  }
  my_fields->grid_dimension[0] = field_size;
  my_fields->grid_end[0] = field_size - 1;

  for (i = 0;i < field_size;i++) {
    my_fields->density[i] = 100.;
    my_fields->HI_density[i] = grackle_data->HydrogenFractionByMass * my_fields->density[i];
    my_fields->HII_density[i] = tiny_number * my_fields->density[i];
    my_fields->HM_density[i] = tiny_number * my_fields->density[i];
    my_fields->HeI_density[i] = (1.0 - grackle_data->HydrogenFractionByMass) *
      my_fields->density[i];
    my_fields->HeII_density[i] = tiny_number * my_fields->density[i];
    my_fields->HeIII_density[i] = tiny_number * my_fields->density[i];
    my_fields->H2I_density[i] = tiny_number * my_fields->density[i];
    my_fields->H2II_density[i] = tiny_number * my_fields->density[i];
    my_fields->DI_density[i] = 2.0 * 3.4e-5 * my_fields->density[i];
    my_fields->DII_density[i] = tiny_number * my_fields->density[i];
    my_fields->HDI_density[i] = tiny_number * my_fields->density[i];
    my_fields->e_density[i] = tiny_number * my_fields->density[i];
    // set metallicity in solar units
    double Zmet = 1.e-4;
    my_fields->metal_density[i] = grackle_data->SolarMetalFractionByMass *
      my_fields->density[i] * Zmet;

    my_fields->x_velocity[i] = 0.0;
    my_fields->y_velocity[i] = 0.0;
    my_fields->z_velocity[i] = 0.0;

    // initilize internal energy (here 1000 K for no reason)
    my_fields->internal_energy[i] = 100. / temperature_units;

    my_fields->volumetric_heating_rate[i] = 0.0;
    my_fields->specific_heating_rate[i] = 1.0;
    my_fields->temperature_floor[i] = 2.73 * (1.+initial_redshift);

    my_fields->isrf_habing[i] = 0.;
    my_fields->RT_HI_ionization_rate[i] = 0.0;
    my_fields->RT_HeI_ionization_rate[i] = 0.0;
    my_fields->RT_HeII_ionization_rate[i] = 0.0;
    my_fields->RT_H2_dissociation_rate[i] = 0.0;
    my_fields->RT_heating_rate[i] = 0.0;
    if (grackle_data->use_dust_evol) {
      double dust2metal = 0.0;
      my_fields->dust_density[i] = dust2metal * my_fields->metal_density[i];
      my_fields->He_gas_metalDensity[i] = (1.-dust2metal) * my_fields->density[i] * grackle_data->SolarAbundances[0] * Zmet;
      my_fields->C_gas_metalDensity[i] = (1.-dust2metal) * my_fields->density[i] * grackle_data->SolarAbundances[1] * Zmet;
      my_fields->N_gas_metalDensity[i] = (1.-dust2metal) * my_fields->density[i] * grackle_data->SolarAbundances[2] * Zmet;
      my_fields->O_gas_metalDensity[i] = (1.-dust2metal) * my_fields->density[i] * grackle_data->SolarAbundances[3] * Zmet;
      my_fields->Ne_gas_metalDensity[i] = (1.-dust2metal) * my_fields->density[i] * grackle_data->SolarAbundances[4] * Zmet;
      my_fields->Mg_gas_metalDensity[i] = (1.-dust2metal) * my_fields->density[i] * grackle_data->SolarAbundances[5] * Zmet;
      my_fields->Si_gas_metalDensity[i] = (1.-dust2metal) * my_fields->density[i] * grackle_data->SolarAbundances[6] * Zmet;
      my_fields->S_gas_metalDensity[i] = (1.-dust2metal) * my_fields->density[i] * grackle_data->SolarAbundances[7] * Zmet;
      my_fields->Ca_gas_metalDensity[i] = (1.-dust2metal) * my_fields->density[i] * grackle_data->SolarAbundances[8] * Zmet;
      my_fields->Fe_gas_metalDensity[i] = (1.-dust2metal) * my_fields->density[i] * grackle_data->SolarAbundances[9] * Zmet;
      my_fields->He_dust_metalDensity[i] = 0.;  // no He in dust
      my_fields->C_dust_metalDensity[i] = dust2metal * my_fields->density[i] * grackle_data->SolarAbundances[1] * Zmet;
      my_fields->N_dust_metalDensity[i] = dust2metal * my_fields->density[i] * grackle_data->SolarAbundances[2] * Zmet;
      my_fields->O_dust_metalDensity[i] = dust2metal * my_fields->density[i] * grackle_data->SolarAbundances[3] * Zmet;
      my_fields->Ne_dust_metalDensity[i] = dust2metal * my_fields->density[i] * grackle_data->SolarAbundances[4] * Zmet;
      my_fields->Mg_dust_metalDensity[i] = dust2metal * my_fields->density[i] * grackle_data->SolarAbundances[5] * Zmet;
      my_fields->Si_dust_metalDensity[i] = dust2metal * my_fields->density[i] * grackle_data->SolarAbundances[6] * Zmet;
      my_fields->S_dust_metalDensity[i] = dust2metal * my_fields->density[i] * grackle_data->SolarAbundances[7] * Zmet;
      my_fields->Ca_dust_metalDensity[i] = dust2metal * my_fields->density[i] * grackle_data->SolarAbundances[8] * Zmet;
      my_fields->Fe_dust_metalDensity[i] = dust2metal * my_fields->density[i] * grackle_data->SolarAbundances[9] * Zmet;
      my_fields->SNe_ThisTimeStep[i] = 1.e-6;  // This is actually the number of SNe per units volume
    }
  }

  return;
}

void copy_fields(grackle_field_data *my_fields, grackle_field_data *old_fields, int field_size ) 
{
  int i;
  double tiny_number = 1.e-20;

  my_fields->grid_dimension[0] = field_size;
  my_fields->grid_start[0] = field_size - 1;
  my_fields->grid_end[0] = field_size - 1;

  for (i = 0;i < field_size;i++) {
    my_fields->density[i] = old_fields->density[i];
    my_fields->HI_density[i] = old_fields->HI_density[i];
    my_fields->HII_density[i] = old_fields->HII_density[i];
    my_fields->HM_density[i] = old_fields->HM_density[i];
    my_fields->HeI_density[i] = old_fields->HeI_density[i];
    my_fields->HeII_density[i] = old_fields->HeII_density[i];
    my_fields->HeIII_density[i] = old_fields->HeIII_density[i];
    my_fields->H2I_density[i] = old_fields->H2I_density[i];
    my_fields->H2II_density[i] = old_fields->H2II_density[i];
    my_fields->DI_density[i] = old_fields->DI_density[i];
    my_fields->DII_density[i] = old_fields->DII_density[i];
    my_fields->HDI_density[i] = old_fields->HDI_density[i];
    my_fields->e_density[i] = old_fields->e_density[i];
    // solar metallicity
    my_fields->metal_density[i] = old_fields->metal_density[i];

    my_fields->x_velocity[i] = 0.0;
    my_fields->y_velocity[i] = 0.0;
    my_fields->z_velocity[i] = 0.0;

    // initilize internal energy (here 1000 K for no reason)
    my_fields->internal_energy[i] = old_fields->internal_energy[i];

    my_fields->volumetric_heating_rate[i] = old_fields->volumetric_heating_rate[i];
    my_fields->specific_heating_rate[i] = old_fields->specific_heating_rate[i];
    my_fields->temperature_floor[i] = old_fields->temperature_floor[i];

    my_fields->isrf_habing[i] = old_fields->isrf_habing[i];
    my_fields->RT_HI_ionization_rate[i] = old_fields->RT_HI_ionization_rate[i];
    my_fields->RT_HeI_ionization_rate[i] = old_fields->RT_HeI_ionization_rate[i];
    my_fields->RT_HeII_ionization_rate[i] = old_fields->RT_HeII_ionization_rate[i];
    my_fields->RT_H2_dissociation_rate[i] = old_fields->RT_H2_dissociation_rate[i];
    my_fields->RT_heating_rate[i] = old_fields->RT_heating_rate[i];

    if (grackle_data->use_dust_evol) {
      double dust2metal = 0.3;
      my_fields->dust_density[i] = old_fields->dust_density[i];
      my_fields->He_gas_metalDensity[i] = old_fields->He_gas_metalDensity[i];
      my_fields->C_gas_metalDensity[i] = old_fields->C_gas_metalDensity[i];
      my_fields->N_gas_metalDensity[i] = old_fields->N_gas_metalDensity[i];
      my_fields->O_gas_metalDensity[i] = old_fields->O_gas_metalDensity[i];
      my_fields->Ne_gas_metalDensity[i] = old_fields->Ne_gas_metalDensity[i];
      my_fields->Mg_gas_metalDensity[i] = old_fields->Mg_gas_metalDensity[i];
      my_fields->Si_gas_metalDensity[i] = old_fields->Si_gas_metalDensity[i];
      my_fields->S_gas_metalDensity[i] = old_fields->S_gas_metalDensity[i];
      my_fields->Ca_gas_metalDensity[i] = old_fields->Ca_gas_metalDensity[i];
      my_fields->Fe_gas_metalDensity[i] = old_fields->Fe_gas_metalDensity[i];
      my_fields->He_dust_metalDensity[i] = old_fields->He_dust_metalDensity[i];
      my_fields->C_dust_metalDensity[i] = old_fields->C_dust_metalDensity[i];
      my_fields->N_dust_metalDensity[i] = old_fields->N_dust_metalDensity[i];
      my_fields->O_dust_metalDensity[i] = old_fields->O_dust_metalDensity[i];
      my_fields->Ne_dust_metalDensity[i] = old_fields->Ne_dust_metalDensity[i];
      my_fields->Mg_dust_metalDensity[i] = old_fields->Mg_dust_metalDensity[i];
      my_fields->Si_dust_metalDensity[i] = old_fields->Si_dust_metalDensity[i];
      my_fields->S_dust_metalDensity[i] = old_fields->S_dust_metalDensity[i];
      my_fields->Ca_dust_metalDensity[i] = old_fields->Ca_dust_metalDensity[i];
      my_fields->Fe_dust_metalDensity[i] = old_fields->Fe_dust_metalDensity[i];
      my_fields->SNe_ThisTimeStep[i] = old_fields->SNe_ThisTimeStep[i];
    }
  }

  return;
}

int main(int argc, char *argv[])
{

  /*********************************************************************
  / Initial setup of units and chemistry objects.
  / This should be done at simulation start.
  *********************************************************************/

  // Enable output
  grackle_verbose = 0;

  // Set initial redshift (for internal units).
  double initial_redshift = 4.;

  // First, set up the units system.
  // These are conversions from code units to cgs.
  code_units my_units;
  my_units.comoving_coordinates = 0; // 1 if cosmological sim, 0 if not
  my_units.density_units = 1.67e-24;
  my_units.length_units = 3.086e21;
  my_units.time_units = 1.0e12;
  my_units.a_units = 1.0; // units for the expansion factor
  // Set expansion factor to 1 for non-cosmological simulation.
  my_units.a_value = 1. / (1. + initial_redshift) / my_units.a_units;
  set_velocity_units(&my_units);

  // set temperature units
  double temperature_units = get_temperature_units(&my_units);

  // Second, create a chemistry object for parameters.  This needs to be a pointer.
  chemistry_data *my_grackle_data;
  my_grackle_data = malloc(sizeof(chemistry_data));
  if (set_default_chemistry_parameters(my_grackle_data) == 0) {
    fprintf(stderr, "Error in set_default_chemistry_parameters.\n");
    return EXIT_FAILURE;
  }

  // Set parameter values for chemistry.
  // Access the parameter storage with the struct you've created
  // or with the grackle_data pointer declared in grackle.h (see further below).
  int i=2;
  if (argc >= 2) sscanf(argv[1],"%d",&i);            
  grackle_data->use_grackle = i; // chemistry on
  grackle_data->with_radiative_cooling = 1; // cooling on
  grackle_data->primordial_chemistry = 2;   // molecular network with H, He, D
  grackle_data->dust_chemistry = 1;         // dust processes
  grackle_data->use_dust_evol = 1;         // track dust formation and destruction
  grackle_data->use_dust_density_field = 1;
  grackle_data->metal_cooling = 1;          // metal cooling on
  grackle_data->UVbackground = 0;           // UV background on
  grackle_data->self_shielding_method = 0;  // self shielding option
  grackle_data->use_specific_heating_rate = 0;  
  grackle_data->accuracy = 0.2;           // fractional accuracy of integration
  grackle_data->max_iterations = 200;           // max iterations
  grackle_data->use_isrf_field = 1;           // flag to turn on ISRF (must provide value in my_fields)
  grackle_data->sne_coeff = 100;           // 
  grackle_data->dust_growth_densref = 2.3e-27;           // 
  grackle_data->grackle_data_file = "../../input/CloudyData_UVB=HM2012.h5"; // data file

  grackle_data->SolarAbundances[0]=0.2485;  // He  (10.93 in units where log[H]=12, so photospheric mass fraction -> Y=0.2485 [Hydrogen X=0.7381]; Anders+Grevesse Y=0.2485, X=0.7314)
  grackle_data->SolarAbundances[1]=2.38e-3; // C   (8.43 -> 2.38e-3, AG=3.18e-3)
  grackle_data->SolarAbundances[2]=0.70e-3; // N   (7.83 -> 0.70e-3, AG=1.15e-3)
  grackle_data->SolarAbundances[3]=5.79e-3; // O   (8.69 -> 5.79e-3, AG=9.97e-3)
  grackle_data->SolarAbundances[4]=1.26e-3; // Ne  (7.93 -> 1.26e-3, AG=1.72e-3)
  grackle_data->SolarAbundances[5]=7.14e-4; // Mg  (7.60 -> 7.14e-4, AG=6.75e-4)
  grackle_data->SolarAbundances[6]=6.71e-3; // Si  (7.51 -> 6.71e-4, AG=7.30e-4)
  grackle_data->SolarAbundances[7]=3.12e-4; // S   (7.12 -> 3.12e-4, AG=3.80e-4)
  grackle_data->SolarAbundances[8]=0.65e-4; // Ca  (6.34 -> 0.65e-4, AG=0.67e-4)
  grackle_data->SolarAbundances[9]=1.31e-3; // Fe (7.50 -> 1.31e-3, AG=1.92e-3)

  // Finally, initialize the chemistry object.
  if (initialize_chemistry_data(&my_units) == 0) {
    fprintf(stderr, "Error in initialize_chemistry_data.\n");
    return EXIT_FAILURE;
  }

  // Create struct for storing grackle field data
  grackle_field_data my_fields;
  int field_size = 1;
  malloc_fields(&my_fields, field_size, grackle_data->use_dust_evol);

  // Initialize values as desired
  double dt = 3.15e7 * 1e6 / my_units.time_units;
  init_fields(&my_fields, grackle_data, my_units, field_size, temperature_units, initial_redshift, dt);

  /*********************************************************************
  / Calling the chemistry solver
  / These routines can now be called during the simulation.
  *********************************************************************/

  // Do the burn-in phase 
  grackle_data->use_grackle = 2; 
  grackle_data->max_iterations = 5;           
  if (solve_chemistry(&my_units, &my_fields, dt) == 0) {
    fprintf(stderr, "Error in solve_chemistry.\n");
    return EXIT_FAILURE;
  }
  fprintf(stderr, "Init Density = %g \n", my_fields.density[0]);
  fprintf(stderr, "Init Electron Density = %g \n", my_fields.e_density[0]);
  fprintf(stderr, "Init HI Density = %g \n", my_fields.HI_density[0]);
  fprintf(stderr, "Init H2 Density = %g %g\n", my_fields.H2I_density[0],my_fields.H2II_density[0]);
  printf("FINISHED BURN-IN FOR %d STEPS USING USE_GRACKLE=%d *************************** \n\n",grackle_data->max_iterations, grackle_data->use_grackle);
  // Store post-burnin data 
  grackle_field_data burnin_fields;
  malloc_fields(&burnin_fields, field_size, grackle_data->use_dust_evol);
  copy_fields(&burnin_fields, &my_fields, field_size);

  // Do Crackle
  fprintf(stderr, "Do crackle...\n");
  grackle_data->use_grackle = 2; 
  grackle_data->max_iterations = 3;
  if (solve_chemistry(&my_units, &my_fields, dt) == 0) {
    fprintf(stderr, "Error in solve_chemistry.\n");
    return EXIT_FAILURE;
  }
  fprintf(stderr, "Density = %g \n", my_fields.density[0]);
  fprintf(stderr, "Electron Density = %g \n", my_fields.e_density[0]);
  fprintf(stderr, "HI Density = %g \n", my_fields.HI_density[0]);
  fprintf(stderr, "H2 Density = %g %g\n", my_fields.H2I_density[0],my_fields.H2II_density[0]);
  if (grackle_data->use_dust_evol) {
    fprintf(stderr, "Dust Density = %g\n", my_fields.dust_density[0]);
  }

  // Calculate temperature.
  gr_float *temperature;
  temperature = malloc(field_size * sizeof(gr_float));
  if (calculate_temperature(&my_units, &my_fields,
                            temperature) == 0) {
    fprintf(stderr, "Error in calculate_temperature.\n");
    return EXIT_FAILURE;
  }
  fprintf(stderr, "Temperature = %g K.\n", temperature[0]);

  //return EXIT_SUCCESS;
  //
  fprintf(stderr, "\nNow do with original grackle...\n");
  grackle_data->max_iterations ++;  // max iterations is offset by 1 in crackle vs grackle
  copy_fields(&my_fields, &burnin_fields, field_size);
  grackle_data->use_grackle = 1; 

  if (solve_chemistry(&my_units, &my_fields, dt) == 0) {
    fprintf(stderr, "Error in solve_chemistry.\n");
    return EXIT_FAILURE;
  }
  fprintf(stderr, "Density = %g \n", my_fields.density[0]);
  fprintf(stderr, "Electron Density = %g \n", my_fields.e_density[0]);
  fprintf(stderr, "HI Density = %g \n", my_fields.HI_density[0]);
  fprintf(stderr, "H2 Density = %g %g\n", my_fields.H2I_density[0],my_fields.H2II_density[0]);
  if (grackle_data->use_dust_evol) fprintf(stderr, "Dust Density = %g\n", my_fields.dust_density[0]);

  // Calculate temperature.
  temperature = malloc(field_size * sizeof(gr_float));
  if (calculate_temperature(&my_units, &my_fields,
                            temperature) == 0) {
    fprintf(stderr, "Error in calculate_temperature.\n");
    return EXIT_FAILURE;
  }

  fprintf(stderr, "Temperature = %g K.\n", temperature[0]);

  // Calculate cooling time.
  gr_float *cooling_time;
  cooling_time = malloc(field_size * sizeof(gr_float));
  if (calculate_cooling_time(&my_units, &my_fields,
                             cooling_time) == 0) {
    fprintf(stderr, "Error in calculate_cooling_time.\n");
    return EXIT_FAILURE;
  }

  fprintf(stderr, "Cooling time = %g s.\n", cooling_time[0] *
          my_units.time_units);

  // Calculate pressure.
  gr_float *pressure;
  double pressure_units = my_units.density_units *
    pow(my_units.velocity_units, 2);
  pressure = malloc(field_size * sizeof(gr_float));
  if (calculate_pressure(&my_units, &my_fields,
                         pressure) == 0) {
    fprintf(stderr, "Error in calculate_pressure.\n");
    return EXIT_FAILURE;
  }

  fprintf(stderr, "Pressure = %le dyne/cm^2.\n", pressure[0]*pressure_units);

  // Calculate gamma.
  gr_float *gamma;
  gamma = malloc(field_size * sizeof(gr_float));
  if (calculate_gamma(&my_units, &my_fields,
                      gamma) == 0) {
    fprintf(stderr, "Error in calculate_gamma.\n");
    return EXIT_FAILURE;
  }

  fprintf(stderr, "gamma = %g.\n", gamma[0]);

  // Calculate dust temperature.
  gr_float *dust_temperature;
  dust_temperature = malloc(field_size * sizeof(gr_float));
  if (calculate_dust_temperature(&my_units, &my_fields,
                      dust_temperature) == 0) {
    fprintf(stderr, "Error in calculate_dust_temperature.\n");
    return EXIT_FAILURE;
  }

  fprintf(stderr, "dust_temperature = %g K.\n", dust_temperature[0]);

  return EXIT_SUCCESS;
}
