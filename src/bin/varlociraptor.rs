// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::process::exit;

use structopt::StructOpt;
use varlociraptor::cli::{Varlociraptor, run};

pub fn main() {
    let opt = Varlociraptor::from_args();

    // setup logger
    fern::Dispatch::new()
        .level(log::LogLevelFilter::Info)
        .chain(std::io::stderr())
        .apply()
        .unwrap();

    exit(match run(opt) {
        Err(e) => {
            println!("Error: {}", e);
            1
        },
        _ => 0
    })
}
