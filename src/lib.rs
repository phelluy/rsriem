/// Physical parameters for the two-fluid isothermal Euler equations.
#[derive(Debug, Clone)]
pub struct EulerIsothermal {
    c: f64,   // speed of sound
    p0: f64,  // reference pressure
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
    let r = (p - prm.p0) / prm.c / prm.c + phi * prm.rz1 + (1. - phi) * prm.rz2;
    [r, r * u, r * v, r * phi]
}

pub fn bal2prim_isot(w: [f64; 4], prm: &EulerIsothermal) -> [f64; 4] {
    let r = w[0];
    let u = w[1] / r;
    let v = w[2] / r;
    let phi = w[3] / r;
    let p = prm.p0 + (r - phi * prm.rz1 - (1. - phi) * prm.rz2) * prm.c * prm.c;
    [p, u, v, phi]
}

/// Physical parameters for the two-fluid Euler equations.
pub struct Euler {
    gamma1: f64, // adiabatic index for fluid 1
    gamma2: f64, // adiabatic index for fluid 2
    pinf1: f64,  // reference pressure for fluid 1
    pinf2: f64,  // reference pressure for fluid 2
}

impl Euler {
    pub fn new(gamma1: f64, gamma2: f64, pinf1: f64, pinf2: f64) -> Self {
        Self {
            gamma1,
            gamma2,
            pinf1,
            pinf2,
        }
    }
    pub fn gamma(&self, phi: f64) -> f64 {
        phi * self.gamma1 + (1. - phi) * self.gamma2
    }
    pub fn pinf(&self, phi: f64) -> f64 {
        phi * self.pinf1 + (1. - phi) * self.pinf2
    }
    pub fn cson(&self, rho: f64, p: f64, phi: f64) -> f64 {
        let gam = self.gamma(phi);
        let pinf = self.pinf(phi);
        (gam * (p + pinf) / rho).sqrt()
    }
}

pub fn prim2bal_euler(y: [f64; 5], prm: &Euler) -> [f64; 5] {
    let r = y[0];
    let u = y[1];
    let v = y[2];
    let p = y[3];
    let phi = y[4];
    let pinf = phi * prm.pinf1 + (1. - phi) * prm.pinf2;
    let gam = phi * prm.gamma1 + (1. - phi) * prm.gamma2;
    let energ = (p + gam * pinf) / (gam - 1.) + 0.5 * r * (u * u + v * v);
    [r, r * u, r * v, energ, r * phi]
}

pub fn bal2prim_euler(w: [f64; 5], prm: &Euler) -> [f64; 5] {
    let r = w[0];
    let u = w[1] / r;
    let v = w[2] / r;
    let energ = w[3];
    let phi = w[4] / r;
    let pinf = phi * prm.pinf1 + (1. - phi) * prm.pinf2;
    let gam = phi * prm.gamma1 + (1. - phi) * prm.gamma2;
    let p = (gam - 1.) * (energ - 0.5 * r * (u * u + v * v)) - gam * pinf;
    [r, u, v, p, phi]
}

/// Physical parameters for the shallow water equations.
#[allow(dead_code)]
pub struct ShallowWater {
    g: f64, // gravity
}

impl ShallowWater {
    pub fn new(g: f64) -> Self {
        Self { g }
    }
}

pub fn prim2bal_sw(y: [f64; 3], _prm: &ShallowWater) -> [f64; 3] {
    let h = y[0];
    let u = y[1];
    let v = y[2];
    [h, h * u, h * v]
}

pub fn bal2prim_sw(w: [f64; 3], _prm: &ShallowWater) -> [f64; 3] {
    let h = w[0];
    let hu = w[1];
    let hv = w[2];
    let u = hu / h;
    let v = hv / h;
    [h, u, v]
}
/// Physical parameters for the Burgers equation.
#[allow(dead_code)]
pub struct Burgers {}

impl Burgers {
    pub fn new() -> Self {
        Self {}
    }
}

pub fn prim2bal_burg(y: [f64; 1], _prm: &Burgers) -> [f64; 1] {
    let u = y[0];
    [u]
}

pub fn bal2prim_burg(w: [f64; 1], _prm: &Burgers) -> [f64; 1] {
    let u = w[0];
    [u]
}

/// Burgers equation
// riemann solver burgers
pub fn riem_burg(wl: [f64; 1], wr: [f64; 1], xi: f64, _prm: &Burgers) -> [f64; 1] {
    let ul = wl[0];
    let ur = wr[0];

    let us = if ul > ur {
        // shock
        let s = 0.5 * (ul + ur);
        if xi < s {
            ul
        } else {
            ur
        }
    } else {
        // rarefaction
        if xi < ul {
            ul
        } else if xi > ur {
            ur
        } else {
            xi
        }
    };
    [us]
}

/// St Venant equations

fn z(a: f64, b: f64, prm: &ShallowWater) -> f64 {
    let g = prm.g;

    let sqrt = f64::sqrt;
    if a > b {
        2. * sqrt(g) / (sqrt(a) + sqrt(b))
    } else {
        sqrt(g) * sqrt((a + b) / 2. / a / b)
    }
}

fn dz(a: f64, b: f64, prm: &ShallowWater) -> f64 {
    let g = prm.g;
    let sqrt = f64::sqrt;
    let pow = f64::powf;
    if a > b {
        -sqrt(g) * pow(sqrt(a) + sqrt(b), -0.2e1) * pow(b, -0.1e1 / 0.2e1)
    } else {
        sqrt(g)
            * sqrt(0.2e1)
            * pow((a + b) / a / b, -0.1e1 / 0.2e1)
            * (0.1e1 / a / b - (a + b) / a * pow(b, -0.2e1))
            / 0.4e1
    }
}

fn k(hl: f64, ul: f64, hr: f64, ur: f64, h: f64, prm: &ShallowWater) -> f64 {
    ul - (h - hl) * z(hl, h, prm) - ur - (h - hr) * z(hr, h, prm)
}

fn dk(hl: f64, hr: f64, h: f64, prm: &ShallowWater) -> f64 {
    -z(hl, h, prm) - (h - hl) * dz(hl, h, prm) - z(hr, h, prm) - (h - hr) * dz(hr, h, prm)
}

// riemann solver st venant
pub fn riem_sw(wl: [f64; 3], wr: [f64; 3], xi: f64, prm: &ShallowWater) -> [f64; 3] {
    let hl = wl[0];
    let ul = wl[1] / hl;
    let vl = wl[2] / hl;
    let hr = wr[0];
    let ur = wr[1] / hr;
    let vr = wr[2] / hr;

    let g = prm.g;

    // critère d'apparition du vide
    let crit = ul - ur + 2. * f64::sqrt(g * hl) + 2. * f64::sqrt(g * hr);

    // méthode de Newton
    let mut hs = 1e-6;
    let mut dh: f64 = 1.;
    let mut iter = 0;

    if crit < 0. {
        hs = 1e-10;
        dh = 0.;
        println!("Vacuum ! crit={}", crit);
    }

    while dh.abs() > 1e-12 {
        dh = -k(hl, ul, hr, ur, hs, prm) / dk(hl, hr, hs, prm);
        hs += dh;
        iter += 1;
        // println!(
        //     "iter ={} hl={} ul={} hr={} ur={} hs={} dh={}",
        //     iter, hl, ul, hr, ur, hs, dh
        // );
        if iter > 20 {
            panic!();
        }
    }
    let sqrt = f64::sqrt;

    let us1 = ul - (hs - hl) * z(hl, hs, prm);
    let us2 = ur + (hs - hr) * z(hr, hs, prm);

    let (lambda1m, lambda1p) = if hs < hl {
        (ul - sqrt(g * hl), us1 - sqrt(g * hs))
    } else {
        let j = -sqrt(g) * sqrt(hl * hs) * sqrt((hl + hs) / 2.);
        (ul + j / hl, ul + j / hl)
    };
    let (lambda2m, lambda2p) = if hs < hr {
        (us2 + sqrt(g * hs), ur + sqrt(g * hr))
    } else {
        let j = sqrt(g) * sqrt(hr * hs) * sqrt((hr + hs) / 2.);
        (ur + j / hr, ur + j / hr)
    };
    // println!("lambda1m={} lambda1p={}", lambda1m, lambda1p);
    // println!("lambda2m={} lambda2p={}", lambda2m, lambda2p);
    // panic!();
    let (h, u, v) = if xi < lambda1m {
        (hl, ul, vl)
    } else if xi < lambda1p {
        let u1 = (ul + 2. * sqrt(g * hl) + 2. * xi) / 3.;
        ((u1 - xi) * (u1 - xi) / g, u1, vl)
    } else if xi < us1 {
        (hs, 0.5 * (us1 + us2), 0.)
    } else if xi < us2 {
        (hs, 0.5 * (us1 + us2), 0.)
    } else if xi < lambda2m {
        (hs, us2, vr)
    } else if xi < lambda2p {
        let u2 = (ur - 2. * sqrt(g * hr) + 2. * xi) / 3.;
        ((u2 - xi) * (u2 - xi) / g, u2, vr)
    } else {
        (hr, ur, vr)
    };
    [h, h * u, h * v]
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

#[allow(clippy::too_many_arguments)]
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
            prm,
        )
        - ul
        - (pl - p)
            * pow(prm.c, -0.2e1)
            * z_isot(
                (p - prm.p0) * pow(prm.c, -0.2e1) + rphi(phil, prm),
                (pl - prm.p0) * pow(prm.c, -0.2e1) + rphi(phil, prm),
                prm,
            )
}

// ignore clippy warning for this func: too many arguments
#[allow(clippy::too_many_arguments)]
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
            prm,
        )
        - (pr - p)
            * pow(cson, -0.4e1)
            * dz_isot(
                (p - prm.p0) * pow(cson, -0.2e1) + rphi(phir, prm),
                (pr - prm.p0) * pow(cson, -0.2e1) + rphi(phir, prm),
                prm,
            )
        + pow(cson, -0.2e1)
            * z_isot(
                (p - prm.p0) * pow(cson, -0.2e1) + rphi(phil, prm),
                (pl - prm.p0) * pow(cson, -0.2e1) + rphi(phil, prm),
                prm,
            )
        - (pl - p)
            * pow(cson, -0.4e1)
            * dz_isot(
                (p - prm.p0) * pow(cson, -0.2e1) + rphi(phil, prm),
                (pl - prm.p0) * pow(cson, -0.2e1) + rphi(phil, prm),
                prm,
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
    let [pl, ul, vl, phil] = yl;
    let yr = bal2prim_isot(wr, prm);
    let [pr, ur, vr, phir] = yr;

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
        //assert!(dp == dp);
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

/// Riemann solver functions for the two-fluid Euler equations.
///

fn phia(pinfa: f64, ga: f64, ta: f64, pa: f64, pl: f64) -> f64 {
    let mut t0;
    t0 = (pl - pa) * (ta - ha(pinfa, ga, ta, pa, pl));
    t0 = t0.abs();
    t0 = t0.sqrt();
    // cas des chocs non entropiques (ne sert qu'au déboguage)
    if pl <= pa {
        t0 = -t0;
    }
    t0
}

fn dphia(pinfa: f64, ga: f64, ta: f64, pa: f64, pl: f64) -> f64 {
    let (t0, pi, pi0);
    pi = pl + pinfa;
    pi0 = pa + pinfa;
    t0 = 0.5 * (2.0 * ta / ((ga - 1.0) * pi0 + (ga + 1.0) * pi)).sqrt()
        - 0.5
            * (((ga - 1.0) * pi0 + (ga + 1.0) * pi) / 2.0 / ta).sqrt()
            * dha(pinfa, ga, ta, pa, pl);
    t0
}

fn ha(pinfa: f64, ga: f64, ta: f64, pa: f64, pl: f64) -> f64 {
    let (t, pi, pi0);
    pi = pl + pinfa;
    pi0 = pa + pinfa;
    // cas du choc
    if pl > pa {
        t = ta * ((ga + 1.0) * pi0 + (ga - 1.0) * pi) / ((ga + 1.0) * pi + (ga - 1.0) * pi0);
    }
    // cas de la détente
    else {
        t = (pi0 / pi).powf(1.0 / ga) * ta;
    }
    t
}

pub fn dha(pinfa: f64, ga: f64, ta: f64, pa: f64, pl: f64) -> f64 {
    let (t, pi, pi0);
    pi = pl + pinfa;
    pi0 = pa + pinfa;
    // choc
    if pl > pa {
        t = -4.0 * ga * ta * pi0 / (pi * (ga + 1.0) + pi0 * (ga - 1.0)).powf(2.0);
    }
    // détente (ne sert qu'au déboguage)
    else {
        t = -ta * (pi0.powf(1.0 / ga)) / ga * (pi.powf(-(ga + 1.0) / ga));
    }
    t
}

fn psia(pinfa: f64, ga: f64, ta: f64, pa: f64, pl: f64) -> f64 {
    let (t1, c0);
    c0 = (ga * (pa + pinfa) * ta).sqrt();
    t1 = 2.0 * c0 / (ga - 1.0) * (((pl + pinfa) / (pa + pinfa)).powf((ga - 1.0) / 2.0 / ga) - 1.0);
    // t1=  2*c0/(ga-1)*(pow((pl+pinfa)/(pa+pinfa),(ga-1.)/2/ga)-1.);
    // println!(
    //     "c0={}, t1={} ga={}, ta={}, pa={}, pl={}",
    //     c0, t1, ga, ta, pa, pl
    // );
    t1
}

fn dpsia(pinfa: f64, ga: f64, ta: f64, pa: f64, pl: f64) -> f64 {
    let (t0, c0);
    c0 = (ga * (pa + pinfa) * ta).sqrt();
    t0 = c0 / ga
        * (pa + pinfa).powf((1.0 - ga) / 2.0 / ga)
        * (pl + pinfa).powf(-(ga + 1.0) / 2.0 / ga);
    //t0=c0/ga*pow(pa+pinfa,(1.-ga)/2./ga)*pow(pl+pinfa,-(ga+1.)/2./ga);

    t0
}

fn pp(pinfa: f64, ga: f64, ta: f64, pa: f64, psi: f64) -> f64 {
    let (mut t0, c, c0);
    c0 = (ga * (pa + pinfa) * ta).sqrt();
    c = (ga - 1.0) / (ga + 1.0) * (psi + 2.0 / (ga - 1.0) * c0);
    t0 = c * c / ga / ta / (pa + pinfa).powf(1.0 / ga);
    t0 = t0.powf(ga / (ga - 1.0)) - pinfa;
    t0
}

fn xhia(pinfa: f64, ga: f64, ta: f64, pa: f64, pl: f64) -> f64 {
    if pl > pa {
        let phia = phia(pinfa, ga, ta, pa, pl);
        //println!("phia={}, ga={}, ta={}, pa={}, pl={}", phia, ga, ta, pa, pl);
        phia
    } else {
        let psia = psia(pinfa, ga, ta, pa, pl);
        //println!("psia={}, ga={}, ta={}, pa={}, pl={}", psia,ga, ta, pa, pl);
        psia
    }
}

fn dxhia(pinfa: f64, ga: f64, ta: f64, pa: f64, pl: f64) -> f64 {
    if pl > pa {
        dphia(pinfa, ga, ta, pa, pl)
    } else {
        dpsia(pinfa, ga, ta, pa, pl)
    }
}

/// Riemann solver for the two-fluid Euler equations.
pub fn riem_euler(wl: [f64; 5], wr: [f64; 5], xi: f64, prm: &Euler) -> [f64; 5] {
    let yl = bal2prim_euler(wl, prm);
    let [rl, ul, vl, pl, phil] = yl;
    let yr = bal2prim_euler(wr, prm);
    let [rr, ur, vr, pr, phir] = yr;

    // println!("yl = {:?}", yl);
    // println!("yr = {:?}", yr);

    let gaml = phil * prm.gamma1 + (1. - phil) * prm.gamma2;
    let gamr = phir * prm.gamma1 + (1. - phir) * prm.gamma2;
    let pinfl = phil * prm.pinf1 + (1. - phil) * prm.pinf2;
    let pinfr = phir * prm.pinf1 + (1. - phir) * prm.pinf2;
    let fl = gaml - 1.;
    let fr = gamr - 1.;

    // println!(
    //     "gaml={} gamr={} pinfl={} pinfr={}",
    //     gaml, gamr, pinfl, pinfr
    // );

    // -p0 is the minimum pressure
    let p0 = f64::min(pinfl, pinfr);

    let mut pn = -p0 + 1e-12;
    let _dvv = xhia(pinfr, fr + 1., 1. / rr, pr, pn) + xhia(pinfl, fl + 1., 1. / rl, pl, pn);
    // println!("delta_v max = {}",dvv/2.);
    // println!("crit_l = {} p0={}",xhia(pinfl,fl+1.,1./rl,pl,pn),p0);
    // println!("crit_r = {} p0={}",xhia(pinfr,fr+1.,1./rr,pr,pn),p0);

    // critère d'apparition du vide
    //let crit = ul-ur-xhia(pinfr,fr+1.,1./rr,pr,-p0+eps)-xhia(pinfl,fl+1.,1./rl,pl,-p0+eps);
    //let crit = ul-ur-xhia(pinfr,fr+1.,1./rr,pr,-p0)-xhia(pinfl,fl+1.,1./rl,pl,-p0);
    let crit =
        (ul - ur - xhia(pinfr, fr + 1., 1. / rr, pr, pn)) - (xhia(pinfl, fl + 1., 1. / rl, pl, pn));
    //println!("crit={}", crit);
    // let crit = 1.;

    // apparition du vide: on prend de la marge
    // à cause des erreurs d'arrondi
    // liés à la fonction xhia (qui contient des fonctions puissance fractionnaire)
    let mut err = f64::MAX;
    if crit < 1e-6 {
        err = 0.;
    }

    let eps = 1e-12;
    let mut iter = 0;
    while err > eps && iter < 100 {
        iter += 1;
        //println!("iter={}", iter);

        let mut ff = ul - ur;
        let mut df = 0.;

        // println!(
        //     "rl={} rr={} ul={} ur={} pl={} pr={} pn={}",
        //     rl, rr, ul, ur, pl, pr, pn
        // );
        ff -= xhia(pinfl, gaml, 1. / rl, pl, pn);
        df -= dxhia(pinfl, gaml, 1. / rl, pl, pn);
        //println!("ff1={} df={}", ff, df);

        // terme de droite (voir Rouy)

        ff -= xhia(pinfr, gamr, 1. / rr, pr, pn);
        df -= dxhia(pinfr, gamr, 1. / rr, pr, pn);
        //println!("ff2={} df={}", ff, df);

        let dp = ff / df;
        //println!("pn={} dp={} err={}", pn, dp, err);

        pn -= dp;
        err = dp.abs(); //pl.max(pr);
    }
    let pm = pn;
    let r2 = 1. / ha(pinfr, gamr, 1. / rr, pr, pn);
    let r1 = 1. / ha(pinfl, gaml, 1. / rl, pl, pn);

    // vitesse de la discontinuité de contact
    // en l'absence de vide, um1 = um2
    let um1 = ul - xhia(pinfl, gaml, 1. / rl, pl, pn);
    let um2 = ur + xhia(pinfr, gamr, 1. / rr, pr, pn);
    let um = if pinfl > pinfr {
        um1
    } else if pinfl < pinfr {
        um2
    } else {
        0.5 * (um1 + um2)
    };

    let mut vit = [0.; 5];
    // vitesses caractéristiques
    if pm <= pl {
        // 1-détente
        vit[0] = ul - 1. / dpsia(pinfl, gaml, 1. / rl, pl, pl) / rl;
        vit[1] = um1 - 1. / dpsia(pinfl, gaml, 1. / rl, pl, pm) / r1;
    } else {
        // choc
        vit[0] =
            ul - f64::sqrt(((gaml + 1.) * (pm + pinfl) + (gaml - 1.) * (pl + pinfl)) * 0.5 / rl);
        vit[1] = vit[0];
    }

    // contact
    vit[2] = um;

    // 3-détente
    if pm <= pr {
        vit[3] = um2 + 1. / dpsia(pinfr, gamr, 1. / rr, pr, pm) / r2;
        vit[4] = ur + 1. / dpsia(pinfr, gamr, 1. / rr, pr, pr) / rr;
    }
    // 3-choc
    else {
        vit[3] =
            ur + f64::sqrt(((gamr + 1.) * (pm + pinfr) + (gamr - 1.) * (pr + pinfr)) * 0.5 / rr);
        vit[4] = vit[3];
    }

    // let mut r;
    // let mut p;
    // // let mut u;
    // let mut f;
    // let mut pinf;
    //let mut phi;
    let (r,u,v,p,phi) =
    // état gauche
    if xi < vit[0] {
        // r = rl;
        // p = pl;
        // u = ul;
        // f = fl;
        // pinf = pinfl;
        //phi = phil;
        (rl, ul, vl,pl,phil)
    } else if xi >= vit[0] && xi < vit[1] {
        let p = pp(pinfl, fl + 1., 1. / rl, pl, ul - xi);
        let r = 1. / ha(pinfl, fl + 1., 1. / rl, pl, p);
        // f = fl;
        // pinf = pinfl;
        let phi = phil;
        let u = ul - psia(pinfl, fl + 1., 1. / rl, pl, p);
        (r,u,vl,p,phi)
    } else if xi >= vit[1] && xi < vit[2] {
        let r = r1;
        let p = pm;   // pm = p0 si vide
        let u = um1;  // um1 = um2 sauf si vide
        // f = fl;
        // pinf = pinfl;
        let phi = phil;
        (r,u,vl,p,phi)
    } else if xi >= vit[2] && xi <= vit[3] {
        let r = r2;
        let p = pm;  // pm = p0 si vide
        let u = um2; // um1 = um2 sauf si vide
        // f = fr;
        // pinf = pinfr;
        let phi = phir;
        (r,u,vr,p,phi)
    } else if xi > vit[3] && xi < vit[4] {
        let p = pp(pinfr, fr + 1., 1. / rr, pr, xi - ur);
        let r = 1. / ha(pinfr, fr + 1., 1. / rr, pr, p);
        // f = fr;
        // pinf = pinfr;
        let phi = phir;
        let u = ur + psia(pinfr, fr + 1., 1. / rr, pr, p);
        (r,u,vr,p,phi)
    } else {
        let r = rr;
        let p = pr;
        let u = ur;
        // f = fr;
        // pinf = pinfr;
        let phi = phir;
        (r,u,vr,p,phi)
    };

    let y = [r, u, v, p, phi];
    prim2bal_euler(y, prm)
}

/// Riemann solver for the two-fluid Euler equations
/// Lagrangian case
pub fn lagflux_euler(
    wl: [f64; 5],
    wr: [f64; 5],
    vn: [f64; 2],
    prm: &Euler,
) -> ([f64; 5], (f64, f64)) {
    let yl = bal2prim_euler(wl, prm);
    let [rl, ul, vl, pl, phil] = yl;
    let yr = bal2prim_euler(wr, prm);
    let [rr, ur, vr, pr, phir] = yr;

    // check that the normal is unitary
    assert!((vn[0] * vn[0] + vn[1] * vn[1] - 1.).abs() < 1e-12);

    let ul = ul * vn[0] + vl * vn[1];
    let ur = ur * vn[0] + vr * vn[1];

    // println!("yl = {:?}", yl);
    // println!("yr = {:?}", yr);

    let gaml = phil * prm.gamma1 + (1. - phil) * prm.gamma2;
    let gamr = phir * prm.gamma1 + (1. - phir) * prm.gamma2;
    let pinfl = phil * prm.pinf1 + (1. - phil) * prm.pinf2;
    let pinfr = phir * prm.pinf1 + (1. - phir) * prm.pinf2;
    let fl = gaml - 1.;
    let fr = gamr - 1.;

    // println!(
    //     "gaml={} gamr={} pinfl={} pinfr={}",
    //     gaml, gamr, pinfl, pinfr
    // );

    // -p0 is the minimum pressure
    let p0 = f64::min(pinfl, pinfr);

    let mut pn = -p0 + 1e-12;
    let _dvv = xhia(pinfr, fr + 1., 1. / rr, pr, pn) + xhia(pinfl, fl + 1., 1. / rl, pl, pn);
    // println!("delta_v max = {}",dvv/2.);
    // println!("crit_l = {} p0={}",xhia(pinfl,fl+1.,1./rl,pl,pn),p0);
    // println!("crit_r = {} p0={}",xhia(pinfr,fr+1.,1./rr,pr,pn),p0);

    // critère d'apparition du vide
    //let crit = ul-ur-xhia(pinfr,fr+1.,1./rr,pr,-p0+eps)-xhia(pinfl,fl+1.,1./rl,pl,-p0+eps);
    //let crit = ul-ur-xhia(pinfr,fr+1.,1./rr,pr,-p0)-xhia(pinfl,fl+1.,1./rl,pl,-p0);
    let crit =
        (ul - ur - xhia(pinfr, fr + 1., 1. / rr, pr, pn)) - (xhia(pinfl, fl + 1., 1. / rl, pl, pn));
    //println!("crit={}", crit);
    // let crit = 1.;

    // apparition du vide: on prend de la marge
    // à cause des erreurs d'arrondi
    // liés à la fonction xhia (qui contient des fonctions puissance fractionnaire)
    let mut err = f64::MAX;
    if crit < 1e-6 {
        err = 0.;
    }

    let eps = 1e-12;
    let mut iter = 0;
    while err > eps && iter < 100 {
        iter += 1;
        //println!("iter={}", iter);

        let mut ff = ul - ur;
        let mut df = 0.;

        // println!(
        //     "rl={} rr={} ul={} ur={} pl={} pr={} pn={}",
        //     rl, rr, ul, ur, pl, pr, pn
        // );
        ff -= xhia(pinfl, gaml, 1. / rl, pl, pn);
        df -= dxhia(pinfl, gaml, 1. / rl, pl, pn);
        //println!("ff1={} df={}", ff, df);

        // terme de droite (voir Rouy)

        ff -= xhia(pinfr, gamr, 1. / rr, pr, pn);
        df -= dxhia(pinfr, gamr, 1. / rr, pr, pn);
        //println!("ff2={} df={}", ff, df);

        let dp = ff / df;
        //println!("pn={} dp={} err={}", pn, dp, err);

        pn -= dp;
        err = dp.abs(); //pl.max(pr);
    }
    // let pm = pn;
    // let r2 = 1. / ha(pinfr, gamr, 1. / rr, pr, pn);
    // let r1 = 1. / ha(pinfl, gaml, 1. / rl, pl, pn);

    // vitesse de la discontinuité de contact
    // en l'absence de vide, um1 = um2
    let um1 = ul - xhia(pinfl, gaml, 1. / rl, pl, pn);
    let um2 = ur + xhia(pinfr, gamr, 1. / rr, pr, pn);
    let um = if pinfl > pinfr {
        um1
    } else if pinfl < pinfr {
        um2
    } else {
        0.5 * (um1 + um2)
    };
    ([0., pn * vn[0], pn * vn[1], pn * um, 0.], (um, pn))
}
/// Riemann solver for the two-fluid Euler equations
/// Lagrangian case
pub fn aleflux_euler(
    wl: [f64; 5],
    wr: [f64; 5],
    vn: [f64; 2],
    prm: &Euler,
) -> ([f64; 5], (f64, f64)) {
    let yl = bal2prim_euler(wl, prm);
    let [rl, ul, vl, pl, phil] = yl;
    let yr = bal2prim_euler(wr, prm);
    let [rr, ur, vr, pr, phir] = yr;

    // check that the normal is unitary
    assert!((vn[0] * vn[0] + vn[1] * vn[1] - 1.).abs() < 1e-12);

    let unl = ul * vn[0] + vl * vn[1];
    let unr = ur * vn[0] + vr * vn[1];

    let utl = -ul * vn[1] + vl * vn[0];
    let utr = -ur * vn[1] + vr * vn[0];

    // println!("yl = {:?}", yl);
    // println!("yr = {:?}", yr);

    let gaml = phil * prm.gamma1 + (1. - phil) * prm.gamma2;
    let gamr = phir * prm.gamma1 + (1. - phir) * prm.gamma2;
    let pinfl = phil * prm.pinf1 + (1. - phil) * prm.pinf2;
    let pinfr = phir * prm.pinf1 + (1. - phir) * prm.pinf2;
    let fl = gaml - 1.;
    let fr = gamr - 1.;

    // println!(
    //     "gaml={} gamr={} pinfl={} pinfr={}",
    //     gaml, gamr, pinfl, pinfr
    // );

    let el = (wl[3] - (wl[1].powi(2) / wl[0] + wl[2].powi(2) / wl[0]) / 2.0) / wl[0];
    let er = (wr[3] - (wr[1].powi(2) / wr[0] + wr[2].powi(2) / wr[0]) / 2.0) / wr[0];
    //let cl = sound_speed(rl, el, phil);
    let cl = (gaml * (pl + pinfl) / rl).sqrt();
    let cr = (gamr * (pr + pinfr) / rr).sqrt();

    let alphal = (gaml + 1.0) / 2.0;
    let alphar = (gamr + 1.0) / 2.0;
    let alpha = alphal.max(alphar);

    let (ar, al);
    if pr > pl {
        let mut a = (pr - pl) / (rr * cr) + unl - unr;
        a = a.max(0.0);
        al = rl * (cl + alpha * a);

        a = (pl - pr) / al + unl - unr;
        a = a.max(0.0);
        ar = rr * (cr + alpha * a);
    } else {
        let mut a = (pl - pr) / (rl * cl) + unl - unr;
        a = a.max(0.0);
        ar = rr * (cr + alpha * a);

        a = (pr - pl) / ar + unl - unr;
        a = a.max(0.0);
        al = rl * (cl + alpha * a);
    }

    let ustar = (al * unl + ar * unr + pl - pr) / (al + ar);
    let pstar = (ar * pl + al * pr - al * ar * (unr - unl)) / (al + ar);
    let r1 = 1.0 / (1.0 / rl + (ar * (unr - unl) + pl - pr) / (al * (al + ar)));
    let r2 = 1.0 / (1.0 / rr + (al * (unr - unl) + pr - pl) / (ar * (al + ar)));
    let e1 = el - pl.powi(2) / (2.0 * al.powi(2)) + pstar.powi(2) / (2.0 * al.powi(2));
    let e2 = er - pr.powi(2) / (2.0 * ar.powi(2)) + pstar.powi(2) / (2.0 * ar.powi(2));

    let u1 = vn[0] * ustar - vn[1] * utl;
    let v1 = vn[1] * ustar + vn[0] * utl;

    let u2 = vn[0] * ustar - vn[1] * utr;
    let v2 = vn[1] * ustar + vn[0] * utr;
    // println!("in u1={} v1={} u2={} v2={}", u1, v1, u2, v2);
    // println!("in r1={} r2={} e1={} e2={}", r1, r2, e1, e2);

    let mut flux = [0.0; 5];
    if (phil - 0.5) * (phir - 0.5) > 0.0 {
        let sigma1 = unl - al / rl;
        let sigma2 = ustar;
        let sigma3 = unr + ar / rr;

        if sigma1 > 0.0 {
            flux[0] = rl * unl;
            flux[1] = rl * ul * unl + pl * vn[0];
            flux[2] = rl * vl * unl + pl * vn[1];
            flux[3] = (wl[3] + pl) * unl;
            flux[4] = rl * unl * phil;
        } else if sigma2 > 0.0 {
            flux[0] = r1 * ustar;
            flux[1] = r1 * u1 * ustar + pstar * vn[0];
            flux[2] = r1 * v1 * ustar + pstar * vn[1];
            flux[3] = (r1 * e1 + r1 * (u1.powi(2) + v1.powi(2)) / 2.0 + pstar) * ustar;
            flux[4] = r1 * ustar * phil;
        } else if sigma3 > 0.0 {
            flux[0] = r2 * ustar;
            flux[1] = r2 * u2 * ustar + pstar * vn[0];
            flux[2] = r2 * v2 * ustar + pstar * vn[1];
            flux[3] = (r2 * e2 + r2 * (u2.powi(2) + v2.powi(2)) / 2.0 + pstar) * ustar;
            flux[4] = r2 * ustar * phir;
        } else {
            flux[0] = rr * unr;
            flux[1] = rr * ur * unr + pr * vn[0];
            flux[2] = rr * vr * unr + pr * vn[1];
            flux[3] = (wr[3] + pr) * unr;
            flux[4] = rr * unr * phir;
        }

        let um = 0.0;
        (flux, (um, pstar))
    } else {
        let um = ustar;
        flux[0] = 0.0;
        flux[1] = pstar * vn[0];
        flux[2] = pstar * vn[1];
        flux[3] = um * pstar;
        flux[4] = 0.0;
        (flux, (ustar, pstar))
    }
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

#[test]
fn test_prim2bal_euler() {
    let gamma1 = 1.4;
    let pinf1 = 0.;
    let gamma2 = 1.4;
    let pinf2 = 0.;
    let prm = Euler::new(gamma1, gamma2, pinf1, pinf2);
    let y = [1., 0., 0., 1., 1.];
    let w = prim2bal_euler(y, &prm);
    println!("w={:?}", w);
    let y2 = bal2prim_euler(w, &prm);
    for i in 0..5 {
        assert!((y[i] - y2[i]).abs() < 1e-12);
    }
    let w2 = prim2bal_euler(y2, &prm);
    println!("w2={:?}", w2);
    for i in 0..5 {
        assert!((w[i] - w2[i]).abs() < 1e-12);
    }
}

#[test]
fn test_lagflux_euler1() {
    let gamma1 = 1.4;
    let pinf1 = 0.;
    let gamma2 = 1.4;
    let pinf2 = 0.;
    let prm = Euler::new(gamma1, gamma2, pinf1, pinf2);
    let w = [1., 1., 0., 3., 1.];
    let y = bal2prim_euler(w, &prm);
    let u = y[1];
    let v = y[2];
    let p = y[3];
    let vn = [1. / (2f64).sqrt(), 1. / (2f64).sqrt()];
    let (flux, (un2, p2)) = lagflux_euler(w, w, vn, &prm);
    println!("flux={:?}", flux);
    let flux2 = [0., p * vn[0], p * vn[1], p * (u * vn[0] + v * vn[1]), 0.];
    for iv in 0..5 {
        assert!((flux[iv] - flux2[iv]).abs() < 1e-12);
    }
    assert!((un2 - u * vn[0] - v * vn[1]).abs() < 1e-12);
    assert!((p2 - p).abs() < 1e-12);
}
#[test]
fn test_lagflux_euler2() {
    let gamma1 = 1.4;
    let pinf1 = 0.;
    let gamma2 = 1.4;
    let pinf2 = 0.;
    let prm = Euler::new(gamma1, gamma2, pinf1, pinf2);
    let w = [1., 1., -1., 3., 1.];
    let y = bal2prim_euler(w, &prm);
    let u = y[1];
    let v = y[2];
    let p = y[3];
    let vn = [1. / (2f64).sqrt(), 1. / (2f64).sqrt()];
    let (flux, (un2, p2)) = lagflux_euler(w, w, vn, &prm);
    println!("flux={:?}", flux);
    let flux2 = [0., p * vn[0], p * vn[1], p * (u * vn[0] + v * vn[1]), 0.];
    for iv in 0..5 {
        assert!((flux[iv] - flux2[iv]).abs() < 1e-12);
    }
    assert!((un2 - u * vn[0] - v * vn[1]).abs() < 1e-12);
    assert!((p2 - p).abs() < 1e-12);
}

#[test]
fn test1_aleflux_euler() {
    let gamma1 = 1.4;
    let pinf1 = 0.;
    let gamma2 = 1.4;
    let pinf2 = 0.;
    let prm = Euler::new(gamma1, gamma2, pinf1, pinf2);
    let w = [1., 1., 1., 3., 1.];
    let y = bal2prim_euler(w, &prm);
    let r = y[0];
    let u = y[1];
    let v = y[2];
    let p = y[3];
    let phi = y[4];
    println!("r={}, u={}, v={}, p={}, phi={}", r, u, v, p, phi);
    let vn = [1. / (2f64).sqrt(), 1. / (2f64).sqrt()];
    //let vn = [1., 0.];
    let (flux, (un3, p2)) = aleflux_euler(w, w, vn, &prm);
    println!("flux={:?}", flux);
    let un2 = u * vn[0] + v * vn[1];
    println!("un2={}, p2={}", un2, p2);
    let flux2 = [
        r * un2,
        r * u * un2 + p * vn[0],
        r * v * un2 + p * vn[1],
        (w[3] + p) * un2,
        r * phi * un2,
    ];
    println!("flux2={:?}", flux2);
    for iv in 0..5 {
        assert!((flux[iv] - flux2[iv]).abs() < 1e-12);
    }
    println!("un2={}, un3={}", un2, un3);
    //assert!((un2 - 0.).abs() < 1e-12);
    assert!((p2 - p).abs() < 1e-12);
}
#[test]
fn test2_aleflux_euler() {
    let gamma1 = 1.4;
    let pinf1 = 0.;
    let gamma2 = 1.4;
    let pinf2 = 0.;
    let prm = Euler::new(gamma1, gamma2, pinf1, pinf2);
    let [rl, ul, vl, pl, phil] = [1., 1., 1., 3., 1.];
    let yl = [rl, ul, vl, pl, phil];
    let [rr, ur, vr, pr, phir] = [1.2, 1., 1., 3., 0.];
    let yr = [rr, ur, vr, pr, phir];
    let wl = prim2bal_euler(yl, &prm);
    let wr = prim2bal_euler(yr, &prm);
    println!("rl={}, ul={}, vl={}, pl={}, phil={}", rl, ul, vl, pl, phil);
    println!("rr={}, ur={}, vr={}, pr={}, phir={}", rr, ur, vr, pr, phir);
    let vn = [1. / (2f64).sqrt(), 1. / (2f64).sqrt()];
    let unl = ul * vn[0] + vl * vn[1];
    //let vn = [1., 0.];
    let (flux, (un2, p2)) = aleflux_euler(wl, wr, vn, &prm);
    println!("flux={:?}", flux);
    let un2 = ul * vn[0] + vl * vn[1];
    println!("un2={}, p2={}", un2, p2);
    let flux2 = [0., pl * vn[0], pl * vn[1], pl * un2, 0.];
    println!("flux2={:?}", flux2);
    for iv in 0..5 {
        assert!((flux[iv] - flux2[iv]).abs() < 1e-12);
    }
    println!("un2={}, unl={}", un2, unl);
    //assert!((un2 - 0.).abs() < 1e-12);
    assert!((p2 - pl).abs() < 1e-12);
}
