/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"

void create_cell_types(void)
{
	// set the random seed
	SeedRandom(parameters.ints("random_seed"));

	/*
	   Put any modifications to default cell definition here if you
	   want to have "inherited" by other cell types.

	   This is a good place to set default functions.
	*/

	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment(&microenvironment);

	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL;
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based;
	cell_defaults.functions.custom_cell_rule = NULL;
	cell_defaults.functions.contact_function = NULL;

	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL;
	cell_defaults.functions.calculate_distance_to_membrane = NULL;

	/*
	   This parses the cell definitions in the XML config file.
	*/

	initialize_cell_definitions_from_pugixml();

	/*
	   Put any modifications to individual cell definitions here.

	   This is a good place to set custom functions.
	*/

	cell_defaults.functions.update_phenotype = phenotype_function;
	cell_defaults.functions.custom_cell_rule = custom_function;
	cell_defaults.functions.contact_function = contact_function;

	/*
	   This builds the map of cell definitions and summarizes the setup.
	*/

	build_cell_definitions_maps();
	display_cell_definitions(std::cout);

	return;
}

void setup_microenvironment(void)
{
	// set domain parameters

	// put any custom code to set non-homogeneous initial conditions or
	// extra Dirichlet nodes here.

	// initialize BioFVM

	initialize_microenvironment();

	return;
}

std::vector<std::vector<double>> create_cell_sphere_positions(double cell_radius, double sphere_radius, std::vector<double> center)
{
	std::vector<std::vector<double>> cells;
	int xc = 0, yc = 0, zc = 0;
	double x_spacing = cell_radius * sqrt(3);
	double y_spacing = cell_radius * 2;
	double z_spacing = cell_radius * sqrt(3);

	std::vector<double> tempPoint(3, 0.0);
	// std::vector<double> cylinder_center(3,0.0);

	double x0 = center[0];
	double y0 = center[1];
	double z0 = center[2];

	for (double z = z0 - sphere_radius; z < z0 + sphere_radius; z += z_spacing, zc++)
	{
		for (double x = x0 - sphere_radius; x < x0 + sphere_radius; x += x_spacing, xc++)
		{
			for (double y = y0 - sphere_radius; y < y0 + sphere_radius; y += y_spacing, yc++)
			{
				tempPoint[0] = x + (zc % 2) * 0.5 * cell_radius;
				tempPoint[1] = y + (xc % 2) * cell_radius;
				tempPoint[2] = z;
				double dr = sqrt((tempPoint[0] - x0) * (tempPoint[0] - x0) + (tempPoint[1] - y0) * (tempPoint[1] - y0) + (tempPoint[2] - z0) * (tempPoint[2] - z0));
				if (dr < sphere_radius)
				{
					cells.push_back(tempPoint);
				}
			}
		}
	}
	return cells;
}

void setup_tissue(void)
{
	// place a cluster of tumor cells at the center

	double cell_radius = cell_defaults.phenotype.geometry.radius;
	double cell_spacing = 0.95 * 2.0 * cell_radius;

	double tumor_radius1 =
		parameters.doubles("first_tumor_radius"); // 250.0;
	double tumor_radius2 =
		parameters.doubles("second_tumor_radius"); // 250.0;

	Cell *pCell = NULL;

	double distance_between_centers = parameters.doubles("distance_between_centers");
	std::cout << "distance between two spheroids is " << distance_between_centers << std::endl;
	
	if (parameters.ints("number_of_spheroids") == 1)
		distance_between_centers = 0.0;

	// first spheroid
	std::vector<double> center1 = {-distance_between_centers * 0.5, 0, 0};
	std::cout << "position of centers1 : " << center1[0] << ' ' << center1[1] << ' ' << center1[2] << std::endl;
	std::vector<std::vector<double>> positions = create_cell_sphere_positions(cell_radius, tumor_radius1, center1);
	std::cout << "creating " << positions.size() << " closely-packed tumor cells in tissue 1... " << std::endl;

	Cell_Definition *pCD = cell_definitions_by_index[0];

	for (int i = 0; i < positions.size(); i++)
	{
		pCell = create_cell(*pCD); // tumor cell
		pCell->assign_position(positions[i]);
	}

	// second spheroid
	if (parameters.ints("number_of_spheroids") == 2)
	{

		std::cout << "# of cell definitions " << cell_definitions_by_index.size() << std::endl;
		Cell_Definition *pCD2;
		if (cell_definitions_by_index.size() > 1)
		{
			pCD2 = cell_definitions_by_index[1];
			std::cout << "use second type of cell for the second spheroid" << std::endl;
		}
		else
		{
			pCD2 = cell_definitions_by_index[0];
			std::cout << "only one cell type exists, so we use the same type for the second spheroid" << std::endl;
		}
		std::vector<double> center2 = {distance_between_centers * 0.5, 0, 0};
		std::cout << "position of centers2 : " << center2[0] << ' ' << center2[1] << ' ' << center2[2] << std::endl;
		positions = create_cell_sphere_positions(cell_radius, tumor_radius2, center2);
		std::cout << "creating " << positions.size() << " closely-packed tumor cells in tissue 2... " << std::endl;

		for (int i = 0; i < positions.size(); i++)
		{
			pCell = create_cell(*pCD2); // tumor cell
			pCell->assign_position(positions[i]);
		}
	}

	// load cells from your CSV file (if enabled)
	// load_cells_from_pugixml();
	return;
}

std::vector<std::string> my_coloring_function(Cell *pCell)
{
	return paint_by_number_cell_coloring(pCell);
}

void phenotype_function(Cell *pCell, Phenotype &phenotype, double dt)
{
	return;
}

void custom_function(Cell *pCell, Phenotype &phenotype, double dt)
{
	return;
}

void contact_function(Cell *pMe, Phenotype &phenoMe, Cell *pOther, Phenotype &phenoOther, double dt)
{
	return;
}