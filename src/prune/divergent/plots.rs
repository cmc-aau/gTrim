use std::collections::HashSet;
use std::fmt::Write as FmtWrite;
use std::path::Path;

use anyhow::{Context, Result};
use resvg::{tiny_skia, usvg};

/// Floor applied to (100 - identity) on log axes, so 100% identity is drawn at 99.9%.
const LOG_FLOOR: f64 = 0.1;

/// Loci with at least this percentage of genomes masked have their panel title highlighted.
const HIGHLIGHT_MASKED_PERC: f64 = 0.5;

/// Scale factor used when rasterizing the SVG to PNG.
const PNG_SCALE: f32 = 1.5;

// figure layout (SVG user units)
const COLS: usize = 7;
const PANEL_W: f64 = 200.0;
const PANEL_H: f64 = 180.0;
const PLOT_LEFT: f64 = 30.0;
const PLOT_RIGHT: f64 = 8.0;
const PLOT_TOP: f64 = 18.0;
const PLOT_BOTTOM: f64 = 18.0;
const MARGIN_LEFT: f64 = 16.0;
const MARGIN_RIGHT: f64 = 8.0;
const HEADER_H: f64 = 46.0;
const FOOTER_H: f64 = 20.0;

const FONT: &str = "DejaVu Sans, Arial, Helvetica, sans-serif";
const GREY: &str = "#9a9a9a";
const RED: &str = "#d62728";
const ORANGE: &str = "#ff7f0e";
const BLUE: &str = "#2f6fbf";

/// Classification of a locus comparison for plotting.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PointClass {
    Normal,
    /// Outlier comparison where the genome locus was not masked by the majority vote.
    Outlier,
    /// Outlier comparison where the genome locus was masked.
    Masked,
}

/// Data required to plot a single locus.
pub struct LocusPlotData {
    pub name: String,
    /// (genome-wide AAI excluding locus, locus identity, class)
    pub points: Vec<(f64, f64, PointClass)>,
    /// (genome-wide AAI excluding locus, expected locus identity), sorted by AAI
    pub trend: Vec<(f64, f64)>,
    pub num_genomes: usize,
    pub num_masked: usize,
    pub num_outliers: usize,
}

/// Axis scale for identity values.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AxisScale {
    Linear,
    /// Log scale on (100 - identity), labelled as identity.
    Log,
}

impl AxisScale {
    fn transform(&self, v: f64) -> f64 {
        match self {
            AxisScale::Linear => v,
            AxisScale::Log => -(100.0 - v).max(LOG_FLOOR).log10(),
        }
    }

    fn ticks(&self, lo: f64) -> Vec<f64> {
        let ticks = match self {
            AxisScale::Linear => {
                let step = if 100.0 - lo > 50.0 { 20.0 } else { 10.0 };
                let mut ticks = Vec::new();
                let mut t = 100.0;
                while t >= lo - 1e-9 {
                    ticks.push(t);
                    t -= step;
                }
                ticks.reverse();
                ticks
            }
            AxisScale::Log => vec![0.0, 50.0, 80.0, 90.0, 95.0, 98.0, 99.0, 99.9],
        };

        ticks.into_iter().filter(|&t| t >= lo - 1e-9).collect()
    }

    fn max_value(&self) -> f64 {
        match self {
            AxisScale::Linear => 100.0,
            AxisScale::Log => 100.0 - LOG_FLOOR,
        }
    }
}

/// Settings shown in the figure header.
pub struct PlotSettings {
    pub title: String,
    pub x_label: String,
    pub y_label: String,
    pub window_size: usize,
    pub outlier_threshold: f64,
    pub total_masked: usize,
    pub total_genome_loci: usize,
}

/// Escape text for inclusion in SVG.
fn escape(text: &str) -> String {
    text.replace('&', "&amp;").replace('<', "&lt;").replace('>', "&gt;")
}

fn format_tick(v: f64) -> String {
    if v.fract().abs() < 1e-9 {
        format!("{}", v as i64)
    } else {
        format!("{}", v)
    }
}

/// Maps identity values to SVG coordinates within a panel's plot area.
struct PanelTransform {
    scale: AxisScale,
    ax: f64,
    ay: f64,
    aw: f64,
    ah: f64,
    tx0: f64,
    tx1: f64,
    ty0: f64,
    ty1: f64,
}

impl PanelTransform {
    fn x(&self, v: f64) -> f64 {
        let t = self.scale.transform(v.min(self.scale.max_value()));
        self.ax + (t - self.tx0) / (self.tx1 - self.tx0) * self.aw
    }

    fn y(&self, v: f64) -> f64 {
        let t = self.scale.transform(v.clamp(0.0, self.scale.max_value()));
        self.ay + self.ah - (t - self.ty0) / (self.ty1 - self.ty0) * self.ah
    }
}

/// Create SVG figure with a scatter plot panel for each locus.
pub fn create_svg(loci: &[LocusPlotData], scale: AxisScale, settings: &PlotSettings) -> String {
    let rows = loci.len().div_ceil(COLS).max(1);
    let width = MARGIN_LEFT + COLS as f64 * PANEL_W + MARGIN_RIGHT;
    let height = HEADER_H + rows as f64 * PANEL_H + FOOTER_H;

    // shared x-axis range across all panels
    let min_x = loci
        .iter()
        .flat_map(|l| l.points.iter().map(|p| p.0))
        .fold(f64::INFINITY, f64::min);
    let x_lo = if min_x.is_finite() { ((min_x / 10.0).floor() * 10.0).max(0.0) } else { 0.0 };
    let y_lo = 0.0;

    let mut svg = String::new();
    let _ = writeln!(
        svg,
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}" font-family="{FONT}">"#,
        w = width,
        h = height
    );
    let _ = writeln!(svg, r#"<rect width="100%" height="100%" fill="white"/>"#);

    // header
    let _ = writeln!(
        svg,
        r#"<text x="{MARGIN_LEFT}" y="18" font-size="13">{}: {} of {} genome loci masked</text>"#,
        escape(&settings.title), settings.total_masked, settings.total_genome_loci
    );
    let mut legend = format!(
        "solid line = running median ({} nearest points), dashed = median \u{2212} {} pts · red = outlier comparisons of masked genome loci, hollow orange = outlier comparisons not masked by majority vote · red title = \u{2265}{}% of genomes masked",
        settings.window_size, settings.outlier_threshold, HIGHLIGHT_MASKED_PERC
    );
    if scale == AxisScale::Log {
        legend.push_str(" · log axes on 100 \u{2212} identity, labelled as identity; 100% drawn at 99.9");
    }
    let _ = writeln!(svg, r##"<text x="{MARGIN_LEFT}" y="34" font-size="9" fill="#444">{}</text>"##, escape(&legend));

    for (idx, locus) in loci.iter().enumerate() {
        let col = idx % COLS;
        let row = idx / COLS;
        let px = MARGIN_LEFT + col as f64 * PANEL_W;
        let py = HEADER_H + row as f64 * PANEL_H;

        let tf = PanelTransform {
            scale,
            ax: px + PLOT_LEFT,
            ay: py + PLOT_TOP,
            aw: PANEL_W - PLOT_LEFT - PLOT_RIGHT,
            ah: PANEL_H - PLOT_TOP - PLOT_BOTTOM,
            tx0: scale.transform(x_lo),
            tx1: scale.transform(scale.max_value()),
            ty0: scale.transform(y_lo),
            ty1: scale.transform(scale.max_value()),
        };

        let x_labels = row == rows - 1 || idx + COLS >= loci.len();
        write_panel(&mut svg, locus, &tf, x_lo, y_lo, settings.outlier_threshold, x_labels, col == 0);
    }

    // axis labels
    let _ = writeln!(
        svg,
        r##"<text x="{}" y="{}" font-size="10" text-anchor="middle" fill="#444">{}</text>"##,
        MARGIN_LEFT + COLS as f64 * PANEL_W / 2.0,
        height - 5.0,
        escape(&settings.x_label)
    );
    let _ = writeln!(
        svg,
        r##"<text transform="translate(10 {}) rotate(-90)" font-size="10" text-anchor="middle" fill="#444">{}</text>"##,
        HEADER_H + rows as f64 * PANEL_H / 2.0,
        escape(&settings.y_label)
    );

    svg.push_str("</svg>\n");
    svg
}

#[allow(clippy::too_many_arguments)]
fn write_panel(
    svg: &mut String,
    locus: &LocusPlotData,
    tf: &PanelTransform,
    x_lo: f64,
    y_lo: f64,
    outlier_threshold: f64,
    x_labels: bool,
    y_labels: bool,
) {
    let (ax, ay, aw, ah) = (tf.ax, tf.ay, tf.aw, tf.ah);

    // title
    let masked_perc = if locus.num_genomes > 0 {
        100.0 * locus.num_masked as f64 / locus.num_genomes as f64
    } else {
        0.0
    };
    let highlight = masked_perc >= HIGHLIGHT_MASKED_PERC;
    let _ = writeln!(
        svg,
        r#"<text x="{}" y="{}" font-size="10" text-anchor="middle"{}>{}</text>"#,
        ax + aw / 2.0,
        ay - 5.0,
        if highlight { format!(r#" fill="{RED}" font-weight="bold""#) } else { String::new() },
        escape(&locus.name)
    );

    // frame, grid lines, and tick labels
    let _ = writeln!(svg, r##"<rect x="{ax:.1}" y="{ay:.1}" width="{aw:.1}" height="{ah:.1}" fill="none" stroke="#cccccc" stroke-width="0.6"/>"##);
    for t in tf.scale.ticks(x_lo) {
        let x = tf.x(t);
        let _ = writeln!(svg, r##"<line x1="{x:.1}" y1="{ay:.1}" x2="{x:.1}" y2="{:.1}" stroke="#ececec" stroke-width="0.5"/>"##, ay + ah);
        if x_labels {
            let _ = writeln!(svg, r##"<text x="{x:.1}" y="{:.1}" font-size="7" text-anchor="middle" fill="#555">{}</text>"##, ay + ah + 9.0, format_tick(t));
        }
    }
    for t in tf.scale.ticks(y_lo) {
        let y = tf.y(t);
        let _ = writeln!(svg, r##"<line x1="{ax:.1}" y1="{y:.1}" x2="{:.1}" y2="{y:.1}" stroke="#ececec" stroke-width="0.5"/>"##, ax + aw);
        if y_labels || tf.scale == AxisScale::Log {
            let _ = writeln!(svg, r##"<text x="{:.1}" y="{:.1}" font-size="7" text-anchor="end" fill="#555">{}</text>"##, ax - 2.0, y + 2.5, format_tick(t));
        }
    }

    // points not identified as outliers, drawing a single point for points that overlap at half-pixel resolution
    let mut seen = HashSet::new();
    let _ = write!(svg, r#"<g fill="{GREY}" fill-opacity="0.45">"#);
    for &(x, y, class) in &locus.points {
        if class != PointClass::Normal {
            continue;
        }
        let (sx, sy) = (tf.x(x), tf.y(y));
        if seen.insert(((sx * 2.0).round() as i64, (sy * 2.0).round() as i64)) {
            let _ = write!(svg, r#"<circle cx="{sx:.1}" cy="{sy:.1}" r="0.9"/>"#);
        }
    }
    svg.push_str("</g>\n");

    // trend and cutoff lines, using a single vertex per pixel column
    let mut trend_path = String::new();
    let mut cutoff_path = String::new();
    let mut last_col = None;
    for (i, &(x, expected)) in locus.trend.iter().enumerate() {
        let sx = tf.x(x);
        let col = sx.round() as i64;
        let is_last = i + 1 == locus.trend.len();
        if Some(col) == last_col && !is_last {
            continue;
        }
        last_col = Some(col);
        let cmd = if trend_path.is_empty() { 'M' } else { 'L' };
        let _ = write!(trend_path, "{cmd}{sx:.1} {:.1} ", tf.y(expected));
        let _ = write!(cutoff_path, "{cmd}{sx:.1} {:.1} ", tf.y(expected - outlier_threshold));
    }
    let _ = writeln!(svg, r#"<path d="{}" fill="none" stroke="{BLUE}" stroke-width="1.1"/>"#, trend_path.trim_end());
    let _ = writeln!(svg, r#"<path d="{}" fill="none" stroke="{BLUE}" stroke-width="0.9" stroke-dasharray="3 2"/>"#, cutoff_path.trim_end());

    // outlier points
    let _ = write!(svg, r#"<g fill="none" stroke="{ORANGE}" stroke-width="0.7">"#);
    for &(x, y, class) in &locus.points {
        if class == PointClass::Outlier {
            let _ = write!(svg, r#"<circle cx="{:.1}" cy="{:.1}" r="1.6"/>"#, tf.x(x), tf.y(y));
        }
    }
    svg.push_str("</g>\n");
    let _ = write!(svg, r#"<g fill="{RED}" fill-opacity="0.85">"#);
    for &(x, y, class) in &locus.points {
        if class == PointClass::Masked {
            let _ = write!(svg, r#"<circle cx="{:.1}" cy="{:.1}" r="1.8"/>"#, tf.x(x), tf.y(y));
        }
    }
    svg.push_str("</g>\n");

    // annotation: bottom-left for linear axes, top-left for log axes
    let ty = match tf.scale {
        AxisScale::Linear => ay + ah - 13.0,
        AxisScale::Log => ay + 9.0,
    };
    let _ = writeln!(
        svg,
        r#"<text x="{:.1}" y="{:.1}" font-size="7"><tspan x="{:.1}">{} masked ({:.2}%)</tspan><tspan x="{:.1}" dy="8">{} outlier comparisons · n {}</tspan></text>"#,
        ax + 3.0,
        ty,
        ax + 3.0,
        locus.num_masked,
        masked_perc,
        ax + 3.0,
        locus.num_outliers,
        locus.points.len()
    );
}

/// Write SVG and PNG figures for the specified axis scale.
pub fn write_plots(
    loci: &[LocusPlotData],
    scale: AxisScale,
    settings: &PlotSettings,
    out_prefix: &Path,
) -> Result<()> {
    let svg = create_svg(loci, scale, settings);

    let svg_file = out_prefix.with_extension("svg");
    std::fs::write(&svg_file, &svg).with_context(|| format!("Failed to create file: {}", svg_file.display()))?;

    let png_file = out_prefix.with_extension("png");
    render_png(&svg, &png_file)?;

    Ok(())
}

/// Rasterize SVG to a PNG file.
fn render_png(svg: &str, png_file: &Path) -> Result<()> {
    let mut opt = usvg::Options::default();
    opt.fontdb_mut().load_system_fonts();

    let tree = usvg::Tree::from_str(svg, &opt).context("Failed to parse generated SVG")?;
    let size = tree
        .size()
        .to_int_size()
        .scale_by(PNG_SCALE)
        .context("Invalid figure size")?;

    let mut pixmap = tiny_skia::Pixmap::new(size.width(), size.height()).context("Failed to allocate PNG image")?;
    resvg::render(&tree, tiny_skia::Transform::from_scale(PNG_SCALE, PNG_SCALE), &mut pixmap.as_mut());
    pixmap
        .save_png(png_file)
        .with_context(|| format!("Failed to create file: {}", png_file.display()))?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_log_transform() {
        let scale = AxisScale::Log;
        assert!((scale.transform(90.0) - -1.0).abs() < 1e-12);
        assert!((scale.transform(99.0) - 0.0).abs() < 1e-12);

        // 100% identity is drawn at the edge of the axis (99.9%)
        assert_eq!(scale.transform(100.0), scale.transform(99.9));
        assert_eq!(scale.transform(100.0), scale.transform(scale.max_value()));
    }

    #[test]
    fn test_linear_ticks() {
        assert_eq!(AxisScale::Linear.ticks(60.0), vec![60.0, 70.0, 80.0, 90.0, 100.0]);
        assert_eq!(AxisScale::Linear.ticks(40.0), vec![40.0, 60.0, 80.0, 100.0]);
        assert_eq!(AxisScale::Linear.ticks(0.0), vec![0.0, 20.0, 40.0, 60.0, 80.0, 100.0]);
        assert_eq!(AxisScale::Log.ticks(60.0), vec![80.0, 90.0, 95.0, 98.0, 99.0, 99.9]);
    }
}
