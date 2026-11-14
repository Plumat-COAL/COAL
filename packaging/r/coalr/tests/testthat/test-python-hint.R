# How coal_pipeline_*() say to get COAL's Python pipeline when it is missing (new-user audit 2026-11-14, U01):
# a development coalr (x.y.z.9000) is not sent to the latest release, whose wheels make codon alignments only.

test_that("a development version names the checkout build, a release the platform wheel", {
  dev <- coalr:::python_install_hint("0.6.0.9000")
  expect_match(dev, "development version", fixed = TRUE)
  expect_match(dev, "COAL_CPP", fixed = TRUE)
  expect_false(grepl("releases/latest", dev, fixed = TRUE))
  rel <- coalr:::python_install_hint("0.6.0")
  expect_match(rel, "https://github.com/plumat-coal/COAL/releases/latest", fixed = TRUE)
  expect_match(rel, "py3-none-<platform>.whl", fixed = TRUE)
  expect_false(grepl("development", rel, fixed = TRUE))
  expect_identical(coalr:::python_install_hint("0.6.1.1"), rel)
  # this package's own version decides by default
  own <- coalr:::python_install_hint()
  expect_identical(own, coalr:::python_install_hint(utils::packageVersion("coalr")))
})
