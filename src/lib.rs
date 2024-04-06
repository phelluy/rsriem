/// Physical parameters for the two-fluid isothermal Euler equations.
pub struct EulerIsothermal {
    c: f64,    // speed of sound
    p0: f64,   // reference pressure
    rz1: f64, // reference density for fluid 1
    rz2: f64, // reference density for fluid 2
}

impl EulerIsothermal {
    pub fn new(c: f64, p0: f64, rz1: f64, rz2: f64) -> Self {
        Self { c, p0, rz1, rz2 }
    }
}

pub fn prim2bal_isot(y: [f64; 4], prm: &EulerIsothermal) -> [f64; 4] {
    let p = y[0];
    let u = y[1];
    let v = y[2];
    let phi = y[3];
    let r = (p-prm.p0)/prm.c/prm.c+phi*prm.rz1+(1.-phi)*prm.rz2;
    [r, r * u, r * v, r * phi]
}

pub fn bal2prim_isot(w: [f64; 4], prm: &EulerIsothermal) -> [f64; 4] {
    let r = w[0];
    let u = w[1]/r;
    let v = w[2]/r;
    let phi = w[3]/r;
    let p = prm.p0 + (r-phi*prm.rz1-(1.-phi)*prm.rz2)*prm.c*prm.c;
    [p, u, v, phi]
}

/// Physical parameters for the two-fluid Euler equations.
struct Euler {
    gamma1: f64, // adiabatic index for fluid 1
    gamma2: f64, // adiabatic index for fluid 2
    pinf1: f64,  // reference pressure for fluid 1
    pinf2: f64,  // reference pressure for fluid 2
}

/// Physical parameters for the shallow water equations.
struct ShallowWater {
    g: f64, // gravity
}

/// Riemann solver functions for the two-fluid isothermal Euler equations.
// Riemann solver
fn z_isot(ra: f64, rb: f64, prm: &EulerIsothermal) -> f64 {
    let ln = f64::ln;
    let pow = f64::powf;
    assert!(ra > 0. && rb > 0.);
    if ra < rb {
        prm.c * (ln(rb) - ln(ra)) / (rb - ra)
    } else {
        prm.c * pow(ra * rb, -0.1e1 / 0.2e1)
    }
}

fn dz_isot(ra: f64, rb: f64, prm: &EulerIsothermal) -> f64 {
    let ln = f64::ln;
    let pow = f64::powf;

    if ra < rb {
        -prm.c / ra / (rb - ra) + prm.c * (ln(rb) - ln(ra)) * pow(rb - ra, -0.2e1)
    } else {
        -prm.c * pow(ra * rb, -0.3e1 / 0.2e1) * rb / 0.2e1
    }
}

fn rphi(phi: f64, prm: &EulerIsothermal) -> f64 {
    phi * prm.rz1 + (1. - phi) * prm.rz2
}

fn pres(r: f64, phi: f64, prm: &EulerIsothermal) -> f64 {
    prm.p0 + (r - rphi(phi, prm)) * prm.c * prm.c
}

fn f_isot(
    pl: f64,
    ul: f64,
    _vl: f64,
    phil: f64,
    pr: f64,
    ur: f64,
    _vr: f64,
    phir: f64,
    p: f64,
    prm: &EulerIsothermal,
) -> f64 {
    let pow = f64::powf;
    ur - (pr - p)
        * pow(prm.c, -0.2e1)
        * z_isot(
            (p - prm.p0) * pow(prm.c, -0.2e1) + rphi(phir, prm),
            (pr - prm.p0) * pow(prm.c, -0.2e1) + rphi(phir, prm),
            prm
        )
        - ul
        - (pl - p)
            * pow(prm.c, -0.2e1)
            * z_isot(
                (p - prm.p0) * pow(prm.c, -0.2e1) + rphi(phil, prm),
                (pl - prm.p0) * pow(prm.c, -0.2e1) + rphi(phil, prm),
                prm
            )
}

fn df_isot(
    pl: f64,
    _ul: f64,
    _vl: f64,
    phil: f64,
    pr: f64,
    _ur: f64,
    _vr: f64,
    phir: f64,
    p: f64,
    prm: &EulerIsothermal,
) -> f64 {
    let pow = f64::powf;
    let cson = prm.c;
    pow(cson, -0.2e1)
        * z_isot(
            (p - prm.p0) * pow(cson, -0.2e1) + rphi(phir, prm),
            (pr - prm.p0) * pow(cson, -0.2e1) + rphi(phir, prm),
            prm
        )
        - (pr - p)
            * pow(cson, -0.4e1)
            * dz_isot(
                (p - prm.p0) * pow(cson, -0.2e1) + rphi(phir, prm),
                (pr - prm.p0) * pow(cson, -0.2e1) + rphi(phir, prm),
                prm
            )
        + pow(cson, -0.2e1)
            * z_isot(
                (p - prm.p0) * pow(cson, -0.2e1) + rphi(phil, prm),
                (pl - prm.p0) * pow(cson, -0.2e1) + rphi(phil, prm),
                prm
            )
        - (pl - p)
            * pow(cson, -0.4e1)
            * dz_isot(
                (p - prm.p0) * pow(cson, -0.2e1) + rphi(phil, prm),
                (pl - prm.p0) * pow(cson, -0.2e1) + rphi(phil, prm),
                prm
            )
}

// riemann solver isothermal model
pub fn riemisot(wl: [f64; 4], wr: [f64; 4], xi: f64, prm: &EulerIsothermal) -> [f64; 4] {
    let rl = wl[0];
    // let ul = wl[1] / rl;
    // let vl = wl[2] / rl;
    // let phil = wl[3] / rl;
    // let pl = pres(rl, phil, prm);
    let rr = wr[0];
    // let ur = wr[1] / rr;
    // let vr = wr[2] / rr;
    // let phir = wr[3] / rr;
    // let pr = pres(rr, phir, prm);
    let yl = bal2prim_isot(wl, prm);
    let [pl, ul, vl, phil]= yl; 
    let yr = bal2prim_isot(wr, prm);
    let [pr, ur, vr, phir]= yr;


    // méthode de Newton
    let mut ps = pres(1e-5, phil, prm).max(pres(1e-5, phir, prm));
    // println!("min pressure={}",ps);
    // let fs = f_isot(pl, ul, vl, phil, pr, ur, vr, phir, ps);
    // println!("f(ps)={}",fs);
    // panic!();
    let mut dp: f64 = 1.;
    let mut iter = 0;
    let iterbound = 40;
    while dp.abs() / prm.p0 > 1e-12 && iter < iterbound {
        let f = f_isot(pl, ul, vl, phil, pr, ur, vr, phir, ps, prm);
        let df = df_isot(pl, ul, vl, phil, pr, ur, vr, phir, ps, prm);
        dp = -f / df;
        ps += dp;
        iter += 1;
        assert!(dp == dp);
        //println!("iter={} dp={} p={} f={} df={}", iter,dp, ps,f,df);
        if iter == iterbound {
            println!("Slow convergence dp={} p={}", dp, ps);
            //panic!();
        }
    }
    let sqrt = f64::sqrt;
    let ra = (ps - prm.p0) / prm.c / prm.c + rphi(phil, prm);
    let rb = (ps - prm.p0) / prm.c / prm.c + rphi(phir, prm);

    let us = ul + (rl - ra) * z_isot(ra, rl, prm);

    let (lambda1m, lambda1p) = if ra < rl {
        (ul - prm.c, us - prm.c)
    } else {
        let lambda = ul - prm.c * sqrt(ra / rl);
        (lambda, lambda)
    };

    let (lambda2m, lambda2p) = if rb < rr {
        (us + prm.c, ur + prm.c)
    } else {
        let lambda = ur + prm.c * sqrt(rb / rr);
        (lambda, us / prm.c * sqrt(rr / rb))
    };

    let (p, u, v, phi) = if xi < lambda1m {
        (pl, ul, vl, phil)
    } else if xi < lambda1p {
        let u1 = xi + prm.c;
        let r1 = rl * ((ul - u1) / prm.c).exp();
        (pres(r1, phil, prm), u1, vl, phil)
    } else if xi < us {
        (ps, us, vl, phil)
    } else if xi < lambda2m {
        (ps, us, vr, phir)
    } else if xi < lambda2p {
        let u2 = xi - prm.c;
        let r2 = rr * ((u2 - ur) / prm.c).exp();
        (pres(r2, phir, prm), u2, vr, phir)
    } else {
        (pr, ur, vr, phir)
    };
    prim2bal_isot([p, u, v, phi], prm)
}

#[cfg(test)]    

// test that bal2prim_isot and prim2bal_isot are consistent
#[test]
fn test_bal2prim_isot() {
    let c = 20.;
    let p0 = 1e5;
    let rz1 = 1.;
    let rz2 = 1000.;
    let prm = EulerIsothermal::new(c, p0, rz1, rz2);
    let y = [1e5, 10., 0., 0.5];
    let w = prim2bal_isot(y, &prm);
    println!("w={:?}", w);
    let y2 = bal2prim_isot(w, &prm);
    for i in 0..4 {
        assert!((y[i] - y2[i]).abs() < 1e-12);
    }
}   

#[test]
fn test_prim2bal_isot() {
    let c = 20.;
    let p0 = 1e5;
    let rz1 = 1.;
    let rz2 = 1000.;
    let prm = EulerIsothermal::new(c, p0, rz1, rz2);
    let w = [1000., 0., 0., 0.];
    let y = bal2prim_isot(w, &prm);
    println!("y={:?}", y);
    let w2 = prim2bal_isot(y, &prm);
    for i in 0..4 {
        assert!((w[i] - w2[i]).abs() < 1e-12);
    }
}  
