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

//! \file
//! \brief Leaky Integrate and Fire neuron type
#ifndef _NEURON_MODEL_LIF_CURR_IMPL_H_
#define _NEURON_MODEL_LIF_CURR_IMPL_H_

#include "neuron_model.h"
#include <common/neuron-typedefs.h>


/////////////////////////////////////////////////////////////
//! definition for LIF neuron parameters
typedef struct neuron_t {
    //! membrane voltage [mV]
    REAL     V_membrane;

    //! membrane resting voltage [mV]
    REAL     V_rest;

    //! membrane resistance [MOhm]
    REAL     R_membrane;

    //! 'fixed' computation parameter - time constant multiplier for
    //! closed-form solution
    //! exp(-(machine time step in ms)/(R * C)) [.]
    REAL     exp_TC;

    //! offset current [nA]
    REAL     I_offset;

    //! countdown to end of next refractory period [timesteps]
    //int32_t  refract_timer;

    //! post-spike reset membrane voltage [mV]
    REAL     V_reset;

    //! refractory time of neuron [timesteps]
    uint32_t  T_refract;

    //! LIF incoming spike times: overwrite after 10 spikes
    float spike_times[10];

    //! LIF incoming spike count
    uint32_t spikeCount;

    //! next internal event time
    float tn;

    //! last event time
    float tl;

    //! Neuron Phases: 0 = Resting,  1 = Subthreshold, 2 = Threshold, 3 = Refractory, 4 = Idle (or Exit), 5 = Error
    uint16_t phase;

    uint16_t waitCounter;

    uint32_t lastProcessedSpikeTime;

    uint32_t lastThresholdTime;

    bool hasSpiked;


} neuron_t;

//! LIF global parameters
typedef struct global_neuron_params_t {
} global_neuron_params_t;


#endif // _NEURON_MODEL_LIF_CURR_IMPL_H_
