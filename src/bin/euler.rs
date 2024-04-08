use rsriem::{bal2prim_euler, prim2bal_euler, riem_euler, Euler};

fn main() {
    // test solver de Riemann Euler bifluide
    let xmin = -10.;
    let xmax = 10.;
    let nx = 1000;
    let tmax = 0.001;

    let rl = 1.;
    let ul = -5.916079783099617;
    //let ul = - 7.;
    //let ul = -4.;
    let vl = 0.;
    let pl = 1.;
    let phil = 1.;

    let rr = 1.;
    //let ur = 4.;
    let ur = 5.916079783099617;
    //let ur = 7.;
    let vr = 0.;
    let pr = 1.;
    let phir = 1.;

    let rl = 1.;
    let ul = -10.;
    //let ul = - 7.;
    //let ul = -4.;
    let vl = 0.;
    let pl = 1e5;
    let phil = 1.;

    let rr = 1000.;
    //let ur = 4.;
    let ur = 10.;
    //let ur = 7.;
    let vr = 0.;
    let pr = 1e5;
    let phir = 0.;


    let gamma1 = 1.4;
    let pinf1 = 0.;
    let gamma2 = 1.4;
    let pinf2 = 0.;

    // params barberon-helluy 2005
    let gamma1 = 1.3;
    let pinf1 = 0.;
    let gamma2 = 3.;
    let pinf2 = 8533e5;
    let prm = Euler::new(gamma1, gamma2, pinf1, pinf2);

    let yl = [rl, ul, vl, pl, phil];
    let yr = [rr, ur, vr, pr, phir];
    let wl = prim2bal_euler(yl, &prm);
    let wr = prim2bal_euler(yr, &prm);

    // solve one riemman problem
    riem_euler(wl, wr, 0., &prm);

    //assert!(1==2);


    let dx = (xmax - xmin) / nx as f64;


    let xi = (0..nx)
        .map(|i| (xmin + (i as f64 + 0.5) * dx))
        .collect::<Vec<f64>>();

    let mut u = vec![0.; nx];
    let mut p = vec![0.; nx];
    let mut rho = vec![0.; nx];
    let mut phi = vec![0.; nx];

    for i in 0..nx {
        let w = riem_euler(wl, wr, xi[i]/(tmax+1e-12), &prm);
        let y = bal2prim_euler(w, &prm);
        u[i] = y[1];
        p[i] = y[3];
        rho[i] = w[0];
        phi[i] = y[4];
    }

    // plot the result
    use rsplot1d::plot;

    println!("rho");
    plot(&xi, &rho, &rho);
    println!("p");
    plot(&xi, &p, &p);
    println!("u");
    plot(&xi, &u, &u);
    println!("phi");
    plot(&xi, &phi, &phi);
}
