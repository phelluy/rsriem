use rsriem::{bal2prim_isot, riemisot, EulerIsothermal, prim2bal_isot};

fn main() {
    // test solver de Riemann isotherme
    let xmin = -10.;
    let xmax = 10.;
    let nx = 1000;
    let tmax = 0.2;

    let pl = 1e5;
    let ul = 0.;
    let vl = 0.;
    let phil = 0.;

    let pr = 1.5e5;
    let ur = 0.;
    let vr = 0.;
    let phir = 1.;

    // let pl = 2.;
    // let ul = 0.;
    // let vl = 0.;
    // let phil = 0.;

    // let pr = 1.;
    // let ur = 0.;
    // let vr = 0.;
    // let phir = 0.;

    let c = 20.;
    let p0 = 1e5;
    let rz1 = 1000.;
    let rz2 = 1000.;
    let prm = EulerIsothermal::new(c, p0, rz1, rz2);

    let yl = [pl, ul, vl, phil];
    let yr = [pr, ur, vr, phir];
    let wl = prim2bal_isot(yl, &prm);
    let wr = prim2bal_isot(yr, &prm);


    let dx = (xmax - xmin) / nx as f64;


    let xi = (0..nx)
        .map(|i| (xmin + (i as f64 + 0.5) * dx))
        .collect::<Vec<f64>>();

    let mut u = vec![0.; nx];
    let mut p = vec![0.; nx];
    let mut rho = vec![0.; nx];

    for i in 0..nx {
        let w = riemisot(wl, wr, xi[i]/(tmax+1e-12), &prm);
        let y = bal2prim_isot(w, &prm);
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
