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

/*! \file
 * \brief implementation of the neuron.h interface.
 */

#include "neuron.h"
#include "neuron_recording.h"
#include "implementations/neuron_impl.h"
#include "plasticity/synapse_dynamics.h"
#include <debug.h>
#include <simulation.h>



//! The key to be used for this core (will be ORed with neuron ID)
static key_t key;

//! A checker that says if this model should be transmitting. If set to false
//! by the data region, then this model should not have a key.
static bool use_key;

//! The number of neurons on the core
static uint32_t n_neurons;

//! The number of clock ticks between sending each spike
static uint32_t time_between_spikes;

//! The expected current clock tick of timer_1 when the next spike can be sent
//static uint32_t expected_time;

//! The recording flags
static uint32_t recording_flags = 0;

static volatile uint32_t exitCounter = 0;

static uint32_t sim_exit_time = 0;

static bool recvd_end_sig = false;
//static bool endsim = false;



//! parameters that reside in the neuron_parameter_data_region
struct neuron_parameters {
    uint32_t timer_start_offset;
    uint32_t time_between_spikes;
    uint32_t has_key;
    uint32_t transmission_key;
    uint32_t n_neurons_to_simulate;
    uint32_t n_synapse_types;
    uint32_t incoming_spike_buffer_size;
};

//! Offset of start of global parameters, in words.
#define START_OF_GLOBAL_PARAMETERS \
    (sizeof(struct neuron_parameters) / sizeof(uint32_t))

//! \brief does the memory copy for the neuron parameters
//! \param[in] address: the address where the neuron parameters are stored
//!     in SDRAM
//! \return bool which is true if the mem copy's worked, false otherwise
static bool neuron_load_neuron_parameters(address_t address) {
    log_debug("loading parameters");
    // call the neuron implementation functions to do the work
    neuron_impl_load_neuron_parameters(
        address, START_OF_GLOBAL_PARAMETERS, n_neurons);
    return true;
}

bool neuron_resume(address_t address) { // EXPORTED
    if (!neuron_recording_reset(n_neurons)){
        log_error("failed to reload the neuron recording parameters");
        return false;
    }

    log_debug("neuron_reloading_neuron_parameters: starting");
    return neuron_load_neuron_parameters(address);
}

bool neuron_initialise(address_t address, address_t recording_address, // EXPORTED
        uint32_t *n_neurons_value, uint32_t *n_synapse_types_value,
        uint32_t *incoming_spike_buffer_size, uint32_t *timer_offset) {
    log_debug("neuron_initialise: starting");
    struct neuron_parameters *params = (void *) address;

    *timer_offset = params->timer_start_offset;
    time_between_spikes = params->time_between_spikes * sv->cpu_clk;
    log_debug("\t back off = %u, time between spikes %u",
            *timer_offset, time_between_spikes);

    // Check if there is a key to use
    use_key = params->has_key;

    // Read the spike key to use
    key = params->transmission_key;

    // output if this model is expecting to transmit
    if (!use_key) {
        //log_info("\tThis model is not expecting to transmit as it has no key");
    } else {
        //log_info("\tThis model is expected to transmit with key = %u", key);
    }

    // Read the neuron details
    n_neurons = params->n_neurons_to_simulate;
    *n_neurons_value = n_neurons;
    *n_synapse_types_value = params->n_synapse_types;

    // Read the size of the incoming spike buffer to use
    *incoming_spike_buffer_size = params->incoming_spike_buffer_size;

    log_info("\t n_neurons = %u, spike buffer size = %u", n_neurons,
            *incoming_spike_buffer_size);

    // Call the neuron implementation initialise function to setup DTCM etc.
    if (!neuron_impl_initialise(n_neurons)) {
        return false;
    }

    // load the data into the allocated DTCM spaces.
    if (!neuron_load_neuron_parameters(address)) {
        return false;
    }

    // setup recording region
    if (!neuron_recording_initialise(
            recording_address, &recording_flags, n_neurons)) {
        return false;
    }

    return true;
}

void neuron_pause(address_t address) { // EXPORTED
    /* Finalise any recordings that are in progress, writing back the final
     * amounts of samples recorded to SDRAM */
    if (recording_flags > 0) {
        log_debug("updating recording regions");
        neuron_recording_finalise();
    }

    // call neuron implementation function to do the work
    neuron_impl_store_neuron_parameters(
            address, START_OF_GLOBAL_PARAMETERS, n_neurons);
}




void neuron_reset_spiked(uint32_t neuron_index){
    neuron_impl_reset_spiked(neuron_index);
}

bool neuron_get_spiked(uint32_t neuron_index){
    return(neuron_impl_spiked(neuron_index));
}

bool neuron_get_lastThresholdTime(uint32_t neuron_index){
    return(neuron_impl_get_lastThresholdTime(neuron_index));
}

uint32_t neuron_update_spiketime(uint32_t time, index_t neuron_index){
    return(neuron_impl_update_spiketime(time, neuron_index));
}

bool neuron_add_spike(uint32_t time, index_t neuron_index){
    if(!neuron_impl_add_spike(neuron_index, time)){
            log_error("Unable to add spike to neuron %u at time = %u", neuron_index, time);
            return false;
    }
    return true;
}



void neuron_set_sim_exit_time(uint32_t time){
    sim_exit_time = time;
    neuron_set_simulation_ticks(time);
}

static void mc_pkt_ignore(uint key, uint payload) {
    use(payload);
    use(payload);
    log_info("Ignore pkt");
  
}

static uint32_t neuron_get_tl(index_t neuron_index){
    return(neuron_impl_get_tl(neuron_index));
}

bool neuron_pdevs_update(){

    input_t external_bias = 0;
    
    for (index_t neuron_index = 0; neuron_index < n_neurons; neuron_index++) {

        //log_info("Doing pdevs update for neuron = %u", neuron_index);
        neuron_recording_setup_for_next_recording();
        neuron_reset_spiked(neuron_index);
    
        bool continueSim = neuron_impl_neuron_update(neuron_index, external_bias,key,use_key);

        if(!continueSim){
            log_info("Calling end_simulation");
            //log_info("Turning off all callbacks at time = %u",  time);
            spin1_callback_off(MCPL_PACKET_RECEIVED);
            //spin1_callback_on(MCPL_PACKET_RECEIVED, mc_pkt_ignore, 1);
            //spin1_callback_off(MC_PACKET_RECEIVED);
            spin1_callback_off(USER_EVENT); 
            spin1_callback_off(DMA_TRANSFER_DONE);
            spin1_delay_us(100);
            uint32_t exitTime = neuron_get_tl(neuron_index);
            neuron_send_terminate_sig(exitTime);
            spin1_delay_us(100);
            
            // simulation_handle_pause_resume(NULL);
            // simulation_ready_to_read();
            //spin1_exit(0);
            //end_sim();
            return false;
        }

    }

    return true;
}

void neuron_set_recvd_end_sig(uint32_t time){
    if((time < sim_exit_time) &&  (time + 15 >= sim_exit_time )){

        recvd_end_sig = true;
        //log_info("Setting recvd_end_sig = %u at time = %u",recvd_end_sig, time);
    }
    
}

void neuron_send_terminate_sig(uint32_t time){

    for (index_t neuron_index = 0; neuron_index < 2; neuron_index++) {
        // call the implementation function (boolean for spike)
        if(use_key){
            log_info("Sending terminate spike (%u, %u), recvd_end_sig= %u", key | neuron_index, time, recvd_end_sig);
            while (!spin1_send_mc_packet(
                            key | neuron_index, time, WITH_PAYLOAD)) {
                        spin1_delay_us(1);
            }
        }
     }       
}

void neuron_add_inputs( // EXPORTED
        index_t synapse_type_index, index_t neuron_index,
        input_t weights_this_timestep) {
    neuron_impl_add_inputs(
            synapse_type_index, neuron_index, weights_this_timestep);
}

#if LOG_LEVEL >= LOG_DEBUG
void neuron_print_inputs(void) { // EXPORTED
    neuron_impl_print_inputs(n_neurons);
}

void neuron_print_synapse_parameters(void) { // EXPORTED
    neuron_impl_print_synapse_parameters(n_neurons);
}
#endif // LOG_LEVEL >= LOG_DEBUG
const char *neuron_get_synapse_type_char(uint32_t synapse_type) { // EXPORTED
    return neuron_impl_get_synapse_type_char(synapse_type);
}

