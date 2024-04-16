use rsriem::{bal2prim_sw, prim2bal_sw, riem_sw, ShallowWater};

fn main() {
    // test solver de Riemann St Venant
    let prm = ShallowWater:: new(9.81);
    let xmin = -10.;
    let xmax = 10.;
    let nx = 1000;
    let _tmax = 1.;

    let hl = 1.;
    let ul = -4.;
    let vl = -1.;

    let hr = 1.;
    let ur = 4.;
    let vr = 1.;

    let yl = [hl, ul, vl];
    let yr = [hr, ur, vr];

    let wl = prim2bal_sw(yl, &prm);
    let wr = prim2bal_sw(yr, &prm);

    // solve one riemman problem
    riem_sw(wl, wr, 0., &prm);

    //assert!(1==2);


    let dx = (xmax - xmin) / nx as f64;


    let xi = (0..nx)
        .map(|i| (xmin + (i as f64 + 0.5) * dx))
        .collect::<Vec<f64>>();

    let mut u = vec![0.; nx];
    let mut v = vec![0.; nx];
    let mut h = vec![0.; nx];
    for i in 0..nx {
        let w = riem_sw(wl, wr, xi[i], &prm);
        let y = bal2prim_sw(w, &prm);
        u[i] = y[1];
        v[i] = y[2];
        h[i] = y[0];

    }

    // plot the result
    use rsplot1d::plot;

    println!("h");
    plot(&xi, &h, &h);
    println!("u");
    plot(&xi, &u, &u);
    println!("v");
    plot(&xi, &v, &v);
}
