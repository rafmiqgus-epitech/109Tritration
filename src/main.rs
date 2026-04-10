use std::env;
use std::fs;
use std::process::ExitCode;

const EXIT_OK: u8 = 0;
const EXIT_ERR: u8 = 84;

#[derive(Clone, Copy, Debug)]
struct Point {
    vol: f64,
    val: f64,
}

fn parse_line(line: &str, lineno: usize) -> Result<Option<Point>, String> {
    let trimmed = line.trim();
    if trimmed.is_empty() {
        return Ok(None);
    }
    let (a, b) = trimmed
        .split_once(';')
        .ok_or_else(|| format!("line {}: expected 'vol;ph'", lineno))?;
    let vol: f64 = a
        .trim()
        .parse()
        .map_err(|_| format!("line {}: invalid volume '{}'", lineno, a.trim()))?;
    let val: f64 = b
        .trim()
        .parse()
        .map_err(|_| format!("line {}: invalid pH '{}'", lineno, b.trim()))?;
    if !vol.is_finite() || !val.is_finite() {
        return Err(format!("line {}: non-finite value", lineno));
    }
    Ok(Some(Point { vol, val }))
}

fn parse_csv(path: &str) -> Result<Vec<Point>, String> {
    let text = fs::read_to_string(path).map_err(|e| format!("{}: {}", path, e))?;
    let mut points = Vec::new();
    for (i, line) in text.lines().enumerate() {
        if let Some(p) = parse_line(line, i + 1)? {
            points.push(p);
        }
    }
    if points.len() < 5 {
        return Err(format!(
            "{}: need at least 5 points, got {}",
            path,
            points.len()
        ));
    }
    for w in points.windows(2) {
        if w[1].vol <= w[0].vol {
            return Err(format!(
                "{}: volumes must be strictly increasing ({} <= {})",
                path, w[1].vol, w[0].vol
            ));
        }
    }
    Ok(points)
}

fn centered_rate(prev: Point, curr: Point, next: Point) -> f64 {
    let h1 = curr.vol - prev.vol;
    let h2 = next.vol - curr.vol;
    let back = (curr.val - prev.val) / h1;
    let fwd = (next.val - curr.val) / h2;
    (h2 * back + h1 * fwd) / (h1 + h2)
}

fn derive_series(points: &[Point]) -> Vec<Point> {
    points
        .windows(3)
        .map(|w| Point {
            vol: w[1].vol,
            val: centered_rate(w[0], w[1], w[2]),
        })
        .collect()
}

fn argmax(series: &[Point]) -> usize {
    let mut best = 0;
    for i in 1..series.len() {
        if series[i].val > series[best].val {
            best = i;
        }
    }
    best
}

fn interpolate_around(series: &[Point], center_idx: usize) -> Result<Vec<Point>, String> {
    if center_idx == 0 || center_idx + 1 >= series.len() {
        return Err("closest point has no neighbours in the second-derivative array".into());
    }
    let prev = series[center_idx - 1];
    let center = series[center_idx];
    let next = series[center_idx + 1];
    let total = ((next.vol - prev.vol) / 0.1).round() as usize;
    let mut out = Vec::with_capacity(total + 1);
    for k in 0..=total {
        let v = prev.vol + (k as f64) * 0.1;
        let y = if v <= center.vol {
            let t = (v - prev.vol) / (center.vol - prev.vol);
            prev.val + t * (center.val - prev.val)
        } else {
            let t = (v - center.vol) / (next.vol - center.vol);
            center.val + t * (next.val - center.val)
        };
        out.push(Point { vol: v, val: y });
    }
    Ok(out)
}

fn first_zero_crossing(samples: &[Point]) -> Result<f64, String> {
    for w in samples.windows(2) {
        if w[0].val > 0.0 && w[1].val <= 0.0 {
            return Ok(w[1].vol);
        }
    }
    Err("no zero crossing found in interpolated second derivative".into())
}

fn round2(x: f64) -> f64 {
    (x * 100.0).round() / 100.0
}

fn print_series(label: &str, series: &[Point]) {
    println!("{}", label);
    for p in series {
        println!("{:.1} ml -> {:.2}", p.vol, round2(p.val));
    }
}

fn print_equivalence(vol: f64) {
    println!("Equivalence point at {:.1} ml", vol);
}

fn print_help() {
    println!("USAGE");
    println!("\t./109titration file");
    println!();
    println!("DESCRIPTION");
    println!("\tfile\ta csv file containing \"vol;ph\" lines");
}

fn run(path: &str) -> Result<(), String> {
    let points = parse_csv(path)?;
    let d1 = derive_series(&points);
    if d1.len() < 3 {
        return Err("not enough points to compute a second derivative".into());
    }
    let closest_vol = d1[argmax(&d1)].vol;
    let d2 = derive_series(&d1);
    let d2_center = match d2.iter().position(|p| (p.vol - closest_vol).abs() < 1e-9) {
        Some(i) => i,
        None => {
            let mut best = 0;
            let mut best_dist = f64::INFINITY;
            for (i, p) in d2.iter().enumerate() {
                let d = (p.vol - closest_vol).abs();
                if d < best_dist {
                    best_dist = d;
                    best = i;
                }
            }
            best
        }
    };
    let interp = interpolate_around(&d2, d2_center)?;
    let final_eq = first_zero_crossing(&interp)?;
    print_series("Derivative:", &d1);
    println!();
    print_equivalence(closest_vol);
    println!();
    print_series("Second derivative:", &d2);
    println!();
    print_series("Second derivative estimated:", &interp);
    println!();
    print_equivalence(final_eq);
    Ok(())
}

fn main() -> ExitCode {
    let args: Vec<String> = env::args().collect();
    if args.len() == 2 && (args[1] == "-h" || args[1] == "--help") {
        print_help();
        return ExitCode::from(EXIT_OK);
    }
    if args.len() != 2 {
        eprintln!("Usage: ./109titration file");
        return ExitCode::from(EXIT_ERR);
    }
    match run(&args[1]) {
        Ok(()) => ExitCode::from(EXIT_OK),
        Err(msg) => {
            eprintln!("{}", msg);
            ExitCode::from(EXIT_ERR)
        }
    }
}
