
//! Logging setup utilities for the application.
//!
//! This module configures logging to both stderr and a log file using the `log4rs` crate.
//! It defines a function to initialize the logger with a consistent format and log level.

use std::path::Path;

use log::LevelFilter;
use log4rs::{
    append::{console::{ConsoleAppender, Target}, file::FileAppender},
    config::{Appender, Config, Root},
    encode::pattern::PatternEncoder,
    filter::threshold::ThresholdFilter,
};

/// Configure logger to write to stderr.
pub fn setup_logger(out_dir: &Path) -> anyhow::Result<()>{
    let level = log::LevelFilter::Info;
    let pattern = "[{d(%Y-%m-%d %H:%M:%S)}] {h({l})}: {m}{n}";

    // log to stderr
    let stderr = ConsoleAppender::builder()
        .encoder(Box::new(PatternEncoder::new(pattern)))
        .target(Target::Stderr)
        .build();

    // log to file
    let logfile = FileAppender::builder()
        .encoder(Box::new(PatternEncoder::new(pattern)))
        .build(out_dir.join("gtrim.log"))?;

    // configure logging
    let config = Config::builder()
        .appender(
            Appender::builder()
                .filter(Box::new(ThresholdFilter::new(level)))
                .build("stderr", Box::new(stderr)),
        ).appender(
            Appender::builder()
                .filter(Box::new(ThresholdFilter::new(level)))
                .build("logfile", Box::new(logfile)),
        )
        .build(Root::builder().appender("stderr").appender("logfile").build(LevelFilter::Trace))
        .expect("Failed to configure logger.");

    log4rs::init_config(config).expect("Failed to initialize logger.");
    Ok(())
}
