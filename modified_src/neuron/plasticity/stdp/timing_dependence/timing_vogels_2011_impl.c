// DO NOT EDIT! THIS FILE WAS GENERATED FROM src/neuron/plasticity/stdp/timing_dependence/timing_vogels_2011_impl.c

/*
 * Copyright (c) 2017-2019 The University of Manchester
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

//! \file
//! \brief Initialisation for timing_vogels_2011_impl.h
#include "timing_vogels_2011_impl.h"

//---------------------------------------
// Globals
//---------------------------------------
//! Lookup table for pre-computed _&tau;_
int16_lut *tau_lookup;

//! Global plasticity parameter data
plasticity_trace_region_data_t plasticity_trace_region_data;

//---------------------------------------
// Functions
//---------------------------------------
address_t timing_initialise(address_t address) {
    log_mini_info("%u", 9207);  /* log_info("timing_initialise: starting");*/
    log_mini_info("%u", 9208);  /* log_info("\tVogels 2011 timing rule");*/

    // Copy parameters
    plasticity_trace_region_data.alpha = (int32_t) address[0];

    // Copy LUTs from following memory
    address_t lut_address = &address[1];
    tau_lookup = maths_copy_int16_lut(&lut_address);

    log_mini_info("%u", 9209);  /* log_info("timing_initialise: completed successfully");*/

    return lut_address;
}
