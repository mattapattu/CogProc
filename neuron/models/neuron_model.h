/*
 * Copyright (c) 2017-2019 The University of Manchester
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/*!
 * \dir
 * \brief Neuronal Soma Models.
 * \file
 * \brief The API for neuron models themselves.
 */

#ifndef _NEURON_MODEL_H_
#define _NEURON_MODEL_H_

#include <common/neuron-typedefs.h>
#include <neuron/threshold_types/threshold_type.h>

//! Forward declaration of neuron type (creates a definition for a pointer to a
//! neuron parameter struct
struct neuron_t;
typedef struct neuron_t neuron_t;
typedef struct neuron_t* neuron_pointer_t;

//! Forward declaration of global neuron parameters
struct global_neuron_params_t;
typedef struct global_neuron_params_t global_neuron_params_t;
//! pointer to global neuron parameters
typedef global_neuron_params_t* global_neuron_params_pointer_t;

//! \brief set the global neuron parameters
//! \param[in] params: The parameters to set
void neuron_model_set_global_neuron_params(
        const global_neuron_params_t *params);

//! \brief primary function called in timer loop after synaptic updates
//! \param[in] num_excitatory_inputs: Number of excitatory receptor types.
//! \param[in] exc_input: Pointer to array of inputs per receptor type received
//!     this timer tick that produce a positive reaction within the neuron in
//!     terms of stimulation.
//! \param[in] num_inhibitory_inputs: Number of inhibitory receptor types.
//! \param[in] inh_input: Pointer to array of inputs per receptor type received
//!     this timer tick that produce a negative reaction within the neuron in
//!     terms of stimulation.
//! \param[in] external_bias: This is the intrinsic plasticity which could be
//!     used for ac, noisy input etc etc. (general purpose input)
//! \param[in,out] neuron: the pointer to a neuron parameter struct which
//!     contains all the parameters for a specific neuron
//! \return the value to be compared with a threshold value to determine if the
//!     neuron has spiked
state_t neuron_model_state_update(int32_t time,
        uint16_t num_excitatory_inputs, const input_t *exc_input,
        uint16_t num_inhibitory_inputs, const input_t *inh_input,
        input_t external_bias, neuron_t *restrict neuron);

//! \brief Indicates that the neuron has spiked
//! \param[in, out] neuron pointer to a neuron parameter struct which contains
//!     all the parameters for a specific neuron
void neuron_model_has_spiked(neuron_t *restrict neuron);

//! \brief get the neuron membrane voltage for a given neuron parameter set
//! \param[in] neuron: a pointer to a neuron parameter struct which contains
//!     all the parameters for a specific neuron
//! \return the membrane voltage for a given neuron with the neuron
//!     parameters specified in neuron
state_t neuron_model_update_membrane_voltage(uint32_t time, neuron_t *neuron);

state_t neuron_model_get_voltage(neuron_t *neuron);



//! \brief printout of state variables i.e. those values that might change
//! \param[in] neuron: a pointer to a neuron parameter struct which contains all
//!     the parameters for a specific neuron
void neuron_model_print_state_variables(const neuron_t *neuron);

//! \brief printout of parameters i.e. those values that don't change
//! \param[in] neuron: a pointer to a neuron parameter struct which contains all
//!     the parameters for a specific neuron
void neuron_model_print_parameters(const neuron_t *neuron);


bool neuron_model_add_spike(neuron_pointer_t neuron, uint32_t  spikeTime);

uint32_t neuron_model_spiketime_pop(neuron_pointer_t neuron);

bool neuron_model_PDevs_sim(neuron_t * neuron, threshold_type_t *threshold_type,  uint32_t nextSpikeTime, key_t key, uint32_t neuron_index, input_t input);

void neuron_model_eit_update(neuron_pointer_t neuron, uint32_t time);


#endif // _NEURON_MODEL_H_
