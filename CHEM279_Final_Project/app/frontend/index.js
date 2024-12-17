const express = require('express'); // Used to handle web requests (post fetch and such)
const path = require('path');       // Path to this index file (our app)
const fs = require('fs').promises;  // File system handler
const WaterSim = require('../backend/build/WaterSim.node'); // Import the C++ Node.js addon
const app = express();              // The web request handler
const port = 3000;                  // The port the server will be attached to.

// Path to the 'err_logs' folder. This is something we should print to with time stamps so that
// if/when we face issues, we can see exactly what they are and when they happened.
const logDir = path.join(__dirname, 'err_logs');

// Ensure the log directory exists
async function ensureLogDir() {
    try {
        await fs.access(logDir);
    } catch (error) {
        await fs.mkdir(logDir, { recursive: true });
    }
}
ensureLogDir();

// Global storage of our points so that they are easier to handle.
let points = [];

app.get('/set/sim/box', async (req, res) => {
    try {
        // WaterSim.SetBoxParams();
        // Send the points as JSON response
        // res.json({ points });
        let w = Number(req.query['w'])
        let h = Number(req.query['h'])
        let d = Number(req.query['d'])
        WaterSim.SetBoxParams(w, h, d)
        res.json({});
    } catch (error) {
        console.error('Error during simulation:', error);
        res.status(500).send('Error executing simulation');
    }
});

// Route for starting the simulation. Accessed as a fetch request from a client.
app.get('/init_simulation', async (req, res) => {
    try {
        let num_part = Number(req.query['num_part']);
        let center = [0.0, 0.0, 0.0];
        let temp = Number(req.query['init_temp']);
        let scale = Number(req.query['sim_scale']);
        let initVel = [0.0, 0.0, 0.0];
        let size = 1.0;
        let jitter_strength = 0.0001;
        const binaryData = WaterSim.SetupSimulation(num_part, scale, center, temp, size, jitter_strength);

        res.set('Content-Type', 'application/octet-stream')
        res.send(binaryData);
    } catch (error) {
        console.error('Error during simulation:', error);
        res.status(500).send('Error executing simulation');
    }
});

// Route for starting the simulation. Accessed as a fetch request from a client.
app.get('/step_sim', async (req, res) => {
    try {
        const binaryData = WaterSim.StepSim();
        res.set('Content-Type', 'application/octet-stream')
        res.send(binaryData);
    } catch (error) {
        console.error('Error during simulation:', error);
        res.status(500).send('Error executing simulation');
    }
});

app.use('/textures/sprites', async (req, res) => {
    res.send(path.join(__dirname, req.url))
});

// Sneds our client the frontend files for rendering.
app.use(express.static(path.join(__dirname, 'public')));

// Starts the server
app.listen(port, () => {
    console.log(`Server running at http://localhost:${port}`);
});
