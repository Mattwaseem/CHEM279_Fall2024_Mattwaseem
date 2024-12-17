// const { response } = require("express");
let curr_positions;
let curr_vels;
let curr_rots;
let curr_angs;
let inv_world_scale = 1.0 / world_scale;

// Initalizes the simulation by calling the C++ script and rendering the generated points.
function initSimulation(num_part, init_temp, sim_scale) {
    // Reset any previous error messages
    document.getElementById("error").textContent = "";

    // Start simulation by fetching data from server
    fetch(`/init_simulation?num_part=${num_part}&init_temp=${init_temp}&sim_scale=${sim_scale}`)
        // Response from server
        .then(response => response.arrayBuffer())
        // Json parse
        .then(data => {
            const view = new DataView(data);
            const num_particles = parameter_container.num_points;
            const pos_start = 0;
            const vel_start = num_particles * 3 * 4;
            const rot_start = vel_start + num_particles * 3 * 4;
            const bond_start = rot_start + num_particles * 4 * 4;

            // This shouldn't happen, but in the event the sim stops sending points over, we will know.
            if (!data || data.byteLength === 0) {
                // console.error('No points data received from server.');
                return;
            }
            // Unpacks the xyz coordinates into an array (best way to stream for a render pass)
            curr_positions = new Float32Array(num_particles * 3);
            curr_rots = new Float32Array(num_particles * 4);
            for (let i = 0; i < num_particles; i++) {

                // C++ Is little endian on ubutnu apparent, so TRUE here is important.
                // I wonder if this can cause an issue if we run the docker on a windows machine
                // instead of mac... mmmmm.
                curr_positions[i * 3] = view.getFloat32(pos_start + i * 12, true) * inv_world_scale;
                curr_positions[i * 3 + 1] = view.getFloat32(pos_start + i * 12 + 4, true) * inv_world_scale;
                curr_positions[i * 3 + 2] = view.getFloat32(pos_start + i * 12 + 8, true) * inv_world_scale;

                curr_rots[i * 4] = view.getFloat32(rot_start + i * 16, true);
                curr_rots[i * 4 + 1] = view.getFloat32(rot_start + i * 16 + 4, true);
                curr_rots[i * 4 + 2] = view.getFloat32(rot_start + i * 16 + 8, true);
                curr_rots[i * 4 + 3] = view.getFloat32(rot_start + i * 16 + 12, true);
            }
            instance_oxy_m.count = num_particles;
            instance_h_m.count = num_particles;

            // Set the position attribute for instancing
            instance_oxy_m.geometry.setAttribute('instancePosition', new THREE.InstancedBufferAttribute(curr_positions, 3));
            instance_h_m.geometry.setAttribute('instancePosition', new THREE.InstancedBufferAttribute(curr_positions, 3));

            instance_oxy_m.geometry.setAttribute('instanceRot', new THREE.InstancedBufferAttribute(curr_rots, 4));
            instance_h_m.geometry.setAttribute('instanceRot', new THREE.InstancedBufferAttribute(curr_rots, 4));
        })
        // Print errors if there are any.
        .catch(error => {
            console.error('Error:', error);
            document.getElementById("error").textContent = "Error loading simulation data";
        });
}


// Step the simulation. This should be renamed from update to step.
function step_sim() {
    fetch('/step_sim')
        // This is the first response gotten back from the websocket of this client js script to the
        // server side node js script
        .then(response => response.arrayBuffer())
        // This happens after the response gets converted into json.
        .then(data => {
            const view = new DataView(data);
            const num_particles = parameter_container.num_points;
            const pos_start = 0;
            const vel_start = num_particles * 3 * 4;
            const rot_start = vel_start + num_particles * 3 * 4;
            const bond_start = rot_start + num_particles * 4 * 4;

            // This shouldn't happen, but in the event the sim stops sending points over, we will know.
            if (!data || data.byteLength === 0) {
                console.error('No points data received from server.');
                return;
            }
            // Unpacks the xyz coordinates into an array (best way to stream for a render pass)
            curr_positions = new Float32Array(num_particles * 3);
            curr_rots = new Float32Array(num_particles * 4);
            for (let i = 0; i < num_particles; i++) {

                // C++ Is little endian on ubutnu apparent, so TRUE here is important.
                // I wonder if this can cause an issue if we run the docker on a windows machine
                // instead of mac... mmmmm.
                curr_positions[i * 3] = view.getFloat32(pos_start + i * 12, true) * inv_world_scale;
                curr_positions[i * 3 + 1] = view.getFloat32(pos_start + i * 12 + 4, true) * inv_world_scale;
                curr_positions[i * 3 + 2] = view.getFloat32(pos_start + i * 12 + 8, true) * inv_world_scale;

                curr_rots[i * 4] = view.getFloat32(rot_start + i * 16, true);
                curr_rots[i * 4 + 1] = view.getFloat32(rot_start + i * 16 + 4, true);
                curr_rots[i * 4 + 2] = view.getFloat32(rot_start + i * 16 + 8, true);
                curr_rots[i * 4 + 3] = view.getFloat32(rot_start + i * 16 + 12, true);
            }
            instance_oxy_m.count = num_particles;
            instance_h_m.count = num_particles;

            // Set the position attribute for instancing
            instance_oxy_m.geometry.setAttribute('instancePosition', new THREE.InstancedBufferAttribute(curr_positions, 3));
            instance_h_m.geometry.setAttribute('instancePosition', new THREE.InstancedBufferAttribute(curr_positions, 3));

            instance_oxy_m.geometry.setAttribute('instanceRot', new THREE.InstancedBufferAttribute(curr_rots, 4));
            instance_h_m.geometry.setAttribute('instanceRot', new THREE.InstancedBufferAttribute(curr_rots, 4));
        })
        // Print errors if there are any.
        .catch(error => {
            console.error('Error:', error);
            document.getElementById("error").textContent = "Error loading simulation data";
        });
}

// Simple boolean flag to let us know when we need to start requesting for simulation steps.
function start_simulation() {
    sim_running = true;
}

function stop_simulation() {
    sim_running = false;
}

// Initialize Sim and threejs
initSimulation(parameter_container.num_points, parameter_container.init_temp);
update_lighting_cols();
update_box_size(parameter_container.num_points);