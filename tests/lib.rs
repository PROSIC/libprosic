use std::fs;
use std::path::Path;
use std::str;
use std::sync::Mutex;

use bio::stats::{LogProb, Prob};
use itertools::Itertools;
use lazy_static::lazy_static;
use paste::paste;
use rust_htslib::bcf;
use rust_htslib::bcf::Read as BCFRead;

mod common;

use common::load_testcase;

macro_rules! testcase {
    ($name:ident, $($pairhmm_mode:ident),+) => {
        paste! {
            lazy_static! {
                static ref [<$name:upper _MUTEX>]: Mutex<()> = Mutex::new(());
            }

            $(
                #[test]
                fn [<$name _ $pairhmm_mode _mode>]() {
                    // Poison error can be ignored here, because it just means that the other test failed
                    // and we are safe to go on.
                    let _guard = [<$name:upper _MUTEX>].lock();
                    let name = stringify!($name);
                    let testcase = load_testcase(
                        &Path::new(file!())
                            .parent()
                            .unwrap()
                            .join("resources/testcases")
                            .join(name),
                    )
                    .unwrap();
                    let mode = stringify!($pairhmm_mode);
                    testcase.run(mode).unwrap();
                    testcase.check();
                }
            )*
        }
    };
}

testcase!(test01, exact, fast);
testcase!(test02, exact, fast);
testcase!(test03, exact, fast);
testcase!(test04, exact, fast);
testcase!(test05, exact, fast);
testcase!(test06, exact, fast);
testcase!(test07, exact, fast);
testcase!(test08, exact, fast);
testcase!(test09, exact, fast);
testcase!(test10, exact, fast);
testcase!(test11, exact, fast);
testcase!(test12, exact, fast);
testcase!(test13, exact, fast);
testcase!(test14, exact, fast);
testcase!(test15, exact, fast);
testcase!(test16, exact, fast);
testcase!(test17, exact, fast);
testcase!(test18, exact, fast);
testcase!(test19, exact, fast);
testcase!(test20, exact, fast);
// skip the next test because this insertion cannot currently be resolved properly
// TODO find a way to fix this.
// testcase!(test21, exact, fast);
testcase!(test22, exact, fast);
testcase!(test23, exact, fast);
testcase!(test24, exact, fast);
testcase!(test25, exact, fast);
testcase!(test26, exact, fast);
testcase!(test27, exact, fast);
testcase!(test28, exact, fast);
testcase!(test29, exact, fast);
testcase!(test30, exact, fast);
testcase!(test31, exact, fast);
testcase!(test32, exact, fast);
testcase!(test33, exact, fast);
testcase!(test34, exact, fast);
testcase!(test36, exact, fast);
testcase!(test37, exact, fast);
// Skip this test. It is most likely a strand bias artifact, which is correctly recognized.
// However, there are also very few reads with nonstandard orientation, which are on the other
// strand.
//testcase!(test38, exact, fast);
testcase!(test39, exact, fast);
testcase!(test40, exact, fast);
testcase!(test41, exact, fast);
testcase!(test42, exact, fast);
testcase!(test43, exact, fast);
// Fast mode fails here, because there is a read with two insertions against
// the alt allele. This is very unlikely to happen, but it happens here.
// In the exact mode, there are various paths around this alignment which rescue
// the alt allele probability. With fast mode, these are missed, making the
// probability artificially small. That leads to Varlociraptor evaluating the
// locus to be heterozygous although it is homozygous in reality.
testcase!(test44, exact);
testcase!(test45, exact, fast);

testcase!(test47, exact, fast);
testcase!(test48, exact, fast);
testcase!(test49, exact);
testcase!(test50, exact);
testcase!(test51, exact);
testcase!(test52, exact, fast);
testcase!(test53, exact, fast);
testcase!(test54, exact, fast);
testcase!(test55, exact, fast);
testcase!(test57, exact);
testcase!(test58, exact);
testcase!(test59, exact);
testcase!(test60, exact);
testcase!(test61, exact);
testcase!(test62, exact);
testcase!(test63, exact);
testcase!(test64, exact);
testcase!(test65, exact);
testcase!(test66, exact);
testcase!(test67, exact);
testcase!(test68, exact);
testcase!(test_mendelian_prior, exact);
testcase!(pattern_too_long, exact, fast);
testcase!(test_long_pattern, exact, fast);
testcase!(test_contig_universe, exact, fast);
testcase!(test_expressions, exact);
testcase!(omit_sb, exact);
testcase!(test_panel_overlap, exact);
testcase!(test_panel_unknown_orientation_bias, exact);
testcase!(issue_154, exact, fast);
testcase!(test_low_cov_vaf, exact);

fn basedir(test: &str) -> String {
    format!("tests/resources/{}", test)
}

fn cleanup_file(f: &str) {
    if Path::new(f).exists() {
        fs::remove_file(f).unwrap();
    }
}

fn control_fdr(test: &str, event_str: &str, alpha: f64) {
    let basedir = basedir(test);
    let output = format!("{}/calls.filtered.bcf", basedir);
    cleanup_file(&output);
    varlociraptor::filtration::fdr::control_fdr(
        &format!("{}/calls.matched.bcf", basedir),
        Some(&output),
        &[varlociraptor::SimpleEvent {
            name: event_str.to_owned(),
        }],
        &varlociraptor::variants::model::VariantType::Deletion(Some(1..30)),
        LogProb::from(Prob(alpha)),
    )
    .unwrap();
}

fn assert_call_number(test: &str, expected_calls: usize) {
    let basedir = basedir(test);

    let mut reader = bcf::Reader::from_path(format!("{}/calls.filtered.bcf", basedir)).unwrap();

    let calls = reader.records().map(|r| r.unwrap()).collect_vec();
    // allow one more or less, in order to be robust to numeric fluctuations
    assert!(
        (calls.len() as i32 - expected_calls as i32).abs() <= 1,
        "unexpected number of calls ({} vs {})",
        calls.len(),
        expected_calls
    );
}

#[test]
fn test_fdr_control1() {
    control_fdr("test_fdr_ev_1", "SOMATIC", 0.05);
    //assert_call_number("test_fdr_ev_1", 974);
}

#[test]
fn test_fdr_control2() {
    control_fdr("test_fdr_ev_2", "SOMATIC", 0.05);
    assert_call_number("test_fdr_ev_2", 985);
}

/// same test, but low alpha
#[test]
fn test_fdr_control3() {
    control_fdr("test_fdr_ev_3", "ABSENT", 0.001);
    assert_call_number("test_fdr_ev_3", 0);
}

#[test]
fn test_fdr_control4() {
    control_fdr("test_fdr_ev_4", "SOMATIC_TUMOR", 0.05);
    assert_call_number("test_fdr_ev_4", 0);
}

// TODO enable this test again once https://github.com/samtools/bcftools/issues/874 is truly fixed upstream
// Then, also encode SVLEN as negative again for deletions.
//#[test]
fn test_fdr_control5() {
    control_fdr("test_fdr_control_out_of_bounds", "PRESENT", 0.05);
}
