// DO NOT EDIT! THIS FILE WAS GENERATED FROM src/neuron/models/neuron_model_lif_impl.c

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
//! \brief Leaky Integrate and Fire neuron implementation
#include "neuron_model_lif_impl.h"
#include <spin1_api.h>
#include <debug.h>


static int32_t last_input_time = 0;
//! \brief simple Leaky I&F ODE
//! \param[in,out] neuron: The neuron to update
//! \param[in] V_prev: previous voltage
//! \param[in] input_this_timestep: The input to apply
static inline void lif_neuron_closed_form(
        neuron_t *neuron, REAL V_prev, input_t input_this_timestep) {
    //REAL alpha = input_this_timestep * neuron->R_membrane + neuron->V_rest;
    // update membrane voltage
    //neuron->V_membrane = alpha - (neuron->exp_TC * (alpha - V_prev));
    int32_t time = spin1_get_simulation_time(); 
    int32_t last_ev_elapsed_time = time - last_input_time;
    //log_info("time = %u, exp(-ms/(RC)) = %11.4k [.]", last_input_time, neuron->exp_TC);
    float decay_factor = 1;
    for(uint32_t k = time ; k > 0; k--){
    	decay_factor = decay_factor*neuron->exp_TC;
    }
    log_mini_info("%u%u%u%x%11.4k", 9130,  time,  last_ev_elapsed_time, float_to_int( decay_factor),  neuron->exp_TC);  /* log_info("time = %u, last_ev_elapsed_time = %u, decay_factor = %f, exp(-ms/(RC)) = %11.4k [.]", time, last_ev_elapsed_time, decay_factor, neuron->exp_TC);*/
    neuron->V_membrane = input_this_timestep * neuron->R_membrane + (neuron->V_membrane/decay_factor);
    log_mini_info("%u%x", 9131, float_to_int( neuron->V_membrane));  /* log_info("V_membrane = %f" , neuron->V_membrane);*/
    last_input_time  = time;
}

void neuron_model_set_global_neuron_params(
        const global_neuron_params_t *params) {
    use(params);
    // Does Nothing - no params
}

state_t neuron_model_state_update(
		uint16_t num_excitatory_inputs, const input_t *exc_input,
		uint16_t num_inhibitory_inputs, const input_t *inh_input,
		input_t external_bias, neuron_t *restrict neuron) {
	//log_info("Exc 1: %12.6k, Exc 2: %12.6k", exc_input[0], exc_input[1]);
	//log_info("Inh 1: %12.6k, Inh 2: %12.6k", inh_input[0], inh_input[1]);

    // If outside of the refractory period
    if (neuron->refract_timer <= 0) {
		REAL total_exc = 0;
		REAL total_inh = 0;

		for (int i=0; i < num_excitatory_inputs; i++) {
			total_exc += exc_input[i];
			//log_info("total_exc: %12.6k, exc_input[i]: %12.6k", total_exc,exc_input[i]);
		}
		for (int i=0; i< num_inhibitory_inputs; i++) {
			total_inh += inh_input[i];
			//log_info("total_inh: %12.6k, inh_input[i]: %12.6k", total_inh,inh_input[i]);
		}
        // Get the input in nA
        input_t input_this_timestep =
                total_exc - total_inh + external_bias + neuron->I_offset;

 log_mini_info("%u%12.6k%12.6k%12.6k", 9132,  input_this_timestep, total_exc, total_inh);  /* log_info("input_this_timestep: %12.6k, total_exc: %12.6k, total_inh: %12.6k", input_this_timestep,total_exc,total_inh);*/

        lif_neuron_closed_form(
                neuron, neuron->V_membrane, input_this_timestep);
    } else {
        // countdown refractory timer
        neuron->refract_timer--;
    }
    return neuron->V_membrane;
}

void neuron_model_has_spiked(neuron_t *restrict neuron) {
    // reset membrane voltage
    neuron->V_membrane = neuron->V_reset;

    // reset refractory timer
    neuron->refract_timer  = neuron->T_refract;
}

state_t neuron_model_get_membrane_voltage(const neuron_t *neuron) {
    return neuron->V_membrane;
}

void neuron_model_print_state_variables(const neuron_t *neuron) {
    log_mini_debug("%u%11.4k", 9133,  neuron->V_membrane);  /* log_debug("V membrane    = %11.4k mv", neuron->V_membrane);*/
}

void neuron_model_print_parameters(const neuron_t *neuron) {
    log_mini_debug("%u%11.4k", 9134,  neuron->V_reset);  /* log_debug("V reset       = %11.4k mv", neuron->V_reset);*/
    log_mini_debug("%u%11.4k", 9135,  neuron->V_rest);  /* log_debug("V rest        = %11.4k mv", neuron->V_rest);*/

    log_mini_debug("%u%11.4k", 9136,  neuron->I_offset);  /* log_debug("I offset      = %11.4k nA", neuron->I_offset);*/
    log_mini_debug("%u%11.4k", 9137,  neuron->R_membrane);  /* log_debug("R membrane    = %11.4k Mohm", neuron->R_membrane);*/

    log_mini_debug("%u%11.4k", 9138,  neuron->exp_TC);  /* log_debug("exp(-ms/(RC)) = %11.4k [.]", neuron->exp_TC);*/

    log_mini_debug("%u%u", 9139,  neuron->T_refract);  /* log_debug("T refract     = %u timesteps", neuron->T_refract);*/
}
