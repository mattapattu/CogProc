// DO NOT EDIT! THIS FILE WAS GENERATED FROM src/neuron/structural_plasticity/synaptogenesis_dynamics_static_impl.c

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

/*!
 * \file
 * \brief This file contains the static (non-)implementation of synaptogenesis.
 * \details No functionality is gained with this class.
 */
#include "synaptogenesis_dynamics.h"
#include <debug.h>

bool synaptogenesis_dynamics_initialise(address_t sdram_sp_address) {
    use(sdram_sp_address);
    return true;
}

bool synaptogenesis_dynamics_rewire(uint32_t time,
        spike_t *spike, address_t *synaptic_row_address, uint32_t *n_bytes) {
    use(time);
    use(spike);
    use(synaptic_row_address);
    use(n_bytes);
    return false;
}

bool synaptogenesis_row_restructure(uint32_t time, address_t row) {
    use(time);
    use(row);
    return false;
}

int32_t synaptogenesis_rewiring_period(void) {
    return -1;
}

bool synaptogenesis_is_fast(void) {
    return false;
}

void synaptogenesis_spike_received(uint32_t time, spike_t spike) {
    use(time);
    use(spike);
}

void print_post_to_pre_entry(void) {
    return;
}
