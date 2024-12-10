use clap::Parser;
use clap::builder::styling::{AnsiColor, Effects, Styles};

fn styles() -> Styles {
    Styles::styled()
        .header(AnsiColor::Yellow.on_default() | Effects::BOLD)
        .usage(AnsiColor::Yellow.on_default() | Effects::BOLD)
        .literal(AnsiColor::Blue.on_default() | Effects::BOLD)
        .placeholder(AnsiColor::Green.on_default())
}
#[derive(Parser, Debug, Clone)]
#[command(version, author, about, long_about = None, styles = styles())]
#[command(
    help_template = "{usage-heading} {usage} \nVersion: {version} {about-section}Author:{author} Email:jiancghen2@genomics.cn/cherryamme@qq.com\n {all-args} {tab}"
)]
pub struct Args {
    /// The path of input file
    #[arg(short, long, num_args = 1..,value_delimiter = ' ')]
    pub inputs: Vec<String>,
    /// The name of outdir
    #[arg(short, long, default_value = "outfile.fq.gz")]
    pub outfile: String,
    /// Nunber of max reads to process
    #[arg(short, long)]
    pub limit: Option<usize>,
    /// Number of threads
    #[arg(short, long, default_value = "10")]
    pub threads: usize,
    /// filter read by min_length
    #[arg(short, long, default_value = "100")]
    pub min_length: usize,
    /// chunk_window
    #[arg(short= 'w', long, default_value = "100")]
    pub chunk_window: usize,
    /// chunk_size
    #[arg(short, long, default_value = "100")]
    pub chunk_size: usize,
    /// correct_ratio
    #[arg(short = 'r', long, default_value = "0.2")]
    pub correct_ratio: f64,
    /// correct_fix_ratio
    #[arg(short = 'f', long, default_value = "0.2")]
    pub correct_fix_ratio: f64,
    /// split_param
    #[arg(long = "id_sep", default_value = "%")]
    pub split_param: char,
    /// type_index
    #[arg(short = 'p', long = "type", default_value = "1")]
    pub type_index: usize,
    /// orient_index
    #[arg(short = 's', long = "strand", default_value = "3")]
    pub orient_index: usize,
}
