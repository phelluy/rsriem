use rsplot1d::plot1d;
use rsriem::{riem_burg, Burgers};

fn main() {
    let xmin = -1.;
    let xmax = 1.;
    let nx = 1000;
    let tmax = 0.5;

    // Riemann problem data: shock wave
    let ul = -1.;
    let ur = 1.;

    let wl = [ul];
    let wr = [ur];

    let prm = Burgers::new();

    let dx = (xmax - xmin) / nx as f64;
    let xi: Vec<f64> = (0..nx).map(|i| xmin + (i as f64 + 0.5) * dx).collect();

    let mut u = vec![0.; nx];

    for i in 0..nx {
        // self-similar solution u(x,t) = W(x/t)
        let w = riem_burg(wl, wr, xi[i] / tmax, &prm);
        u[i] = w[0];
    }

    // exact solution plot
    // plot both as the same because it is exact
    plot1d(&xi, &u, &u);
}
