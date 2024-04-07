use rsriem::{bal2prim_euler, prim2bal_euler, riem_euler, Euler};

fn main() {
    // test solver de Riemann Euler bifluide
    let xmin = -1.;
    let xmax = 1.;
    let nx = 1000;
    let tmax = 0.2;

    let rl = 1.;
    let ul = -0.5;
    let vl = 0.;
    let pl = 1.;
    let phil = 1.;

    let rr = 1.;
    let ur = 0.5;
    let vr = 0.;
    let pr = 1.;
    let phir = 1.;


    let gamma1 = 1.4;
    let pinf1 = 0.;
    let gamma2 = 1.4;
    let pinf2 = 0.;
    let prm = Euler::new(gamma1, gamma2, pinf1, pinf2);

    let yl = [rl, ul, vl, pl, phil];
    let yr = [rr, ur, vr, pr, phir];
    let wl = prim2bal_euler(yl, &prm);
    let wr = prim2bal_euler(yr, &prm);

    // solve one riemman problem
    riem_euler(wl, wr, 0., &prm);



    let dx = (xmax - xmin) / nx as f64;


    let xi = (0..nx)
        .map(|i| (xmin + (i as f64 + 0.5) * dx))
        .collect::<Vec<f64>>();

    let mut u = vec![0.; nx];
    let mut p = vec![0.; nx];
    let mut rho = vec![0.; nx];

    for i in 0..nx {
        let w = riem_euler(wl, wr, xi[i]/(tmax+1e-12), &prm);
        let y = bal2prim_euler(w, &prm);
        u[i] = y[1];
        p[i] = y[0];
        rho[i] = w[0];
    }

    // plot the result
    use rsplot1d::plot;

    plot(&xi, &rho, &rho);
    plot(&xi, &p, &p);
    plot(&xi, &u, &u);
}
