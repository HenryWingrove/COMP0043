extern crate rand_distr;

struct GeometricBrownianMotion {
    spot: f64,
    vol: f64,
    r: f64,
    q: f64,
    dt: f64
}

impl GeometricBrownianMotion {
    fn get_sde(&self, nsteps: u32) {
        let mut dx = vec![1.0; nsteps];
        let dist = rand_distr::Normal::new(0.0, 1.0);
        dx[0] = self.spot;
        for i in 0..nsteps {
            let w = rng.sample(dist);
            let diffusion = (self.vol.pow(i) * W);
            let x = self.dt * (self.r - 0.5*self.vol.powi(2)) * diffusion.sqrt();
            dx[i] = dx[i-1] * x.exp()
        }
    }
}