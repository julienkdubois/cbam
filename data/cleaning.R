suppressPackageStartupMessages({
	library(dplyr)
	library(readr)
	library(tidyr)
	library(stringr)
	library(lubridate)
	library(purrr)
	library(countrycode)
})

library(ecb)
# Taux directeur BCE deposit facility
ecb_rate <- get_data("FM.Q.U2.EUR.RT.MM.EURIBOR3MD_.HSTA")
# ou mieux : taux de dépôt
ecb_rate <- get_data("FM.B.U2.EUR.4F.KR.DFR.LEV")

# Return the first matching column found in a data.frame.
pick_first_col <- function(df, candidates) {
	col <- intersect(candidates, names(df))
	if (length(col) == 0) return(rep(NA, nrow(df)))
	df[[col[1]]]
}

# Prepare ECB policy rate as quarterly series for EU27 fallback.
clean_ecb_policy_rate <- function(df) {
	if (is.null(df) || nrow(df) == 0) {
		return(tibble(date = character(), ecb_policy_rate = numeric()))
	}

	ecb_raw <- tibble(
		date_change = as.Date(pick_first_col(df, c("obstime", "TIME_PERIOD", "TIME", "time", "period", "time_period", "TIMEPERIOD"))),
		rate = suppressWarnings(as.numeric(pick_first_col(df, c("obsvalue", "OBS_VALUE", "obs_value", "value", "OBS"))))
	) %>%
		filter(!is.na(date_change), is.finite(rate)) %>%
		arrange(date_change)

	if (nrow(ecb_raw) == 0) {
		return(tibble(date = character(), ecb_policy_rate = numeric()))
	}

	start_q <- floor_date(min(ecb_raw$date_change), "quarter")
	end_q <- floor_date(max(ecb_raw$date_change), "quarter")

	quarters_grid <- tibble(
		date_end = seq.Date(start_q, end_q, by = "quarter") + months(3) - days(1)
	) %>%
		mutate(date = paste0(year(date_end), "Q", quarter(date_end)))

	idx <- findInterval(quarters_grid$date_end, ecb_raw$date_change)

	quarters_grid %>%
		mutate(ecb_policy_rate = if_else(idx > 0, ecb_raw$rate[idx], NA_real_)) %>%
		select(date, ecb_policy_rate)
}

# Resolve paths from this script location so execution does not depend on getwd().
get_script_dir <- function() {
	args <- commandArgs(trailingOnly = FALSE)
	file_arg <- "--file="
	idx <- grep(file_arg, args)

	if (length(idx) > 0) {
		script_path <- normalizePath(sub(file_arg, "", args[idx[1]]), winslash = "/", mustWork = FALSE)
		return(dirname(script_path))
	}

	if (!is.null(sys.frames()[[1]]$ofile)) {
		return(dirname(normalizePath(sys.frames()[[1]]$ofile, winslash = "/", mustWork = FALSE)))
	}

	normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

script_dir <- get_script_dir()

# Resolve the project data root by testing common locations.
resolve_data_root <- function() {
	candidates <- c(
		script_dir,
		file.path(script_dir, "data"),
		file.path(getwd(), "cbam", "data"),
		file.path(getwd(), "data")
	)

	for (p in candidates) {
		if (dir.exists(file.path(p, "raw"))) {
			return(normalizePath(p, winslash = "/", mustWork = FALSE))
		}
	}

	# Conservative fallback requested by user: <cwd>/cbam/data
	normalizePath(file.path(getwd(), "cbam", "data"), winslash = "/", mustWork = FALSE)
}

data_root <- resolve_data_root()
raw_dir <- file.path(data_root, "raw")
out_dir <- file.path(data_root, "clean")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Keep all countries/areas available in the sources; do not hard-filter to a fixed list.

# Convert country labels/codes to ISO3.
to_iso3 <- function(x) {
	x <- str_trim(as.character(x))

	custom <- c(
		"Korea" = "KOR",
		"Korea, Rep." = "KOR",
		"Turkiye" = "TUR",
		"Turkey" = "TUR",
		"Czechia" = "CZE",
		"Czech Republic" = "CZE",
		"Slovak Republic" = "SVK",
		"United States" = "USA",
		"United Kingdom" = "GBR",
		"European Union" = "EU27",
		"EU" = "EU27",
		"EU27_2020" = "EU27"
	)

	out <- rep(NA_character_, length(x))

	# Step 1: Detect EU/Euro Area variants (any pattern containing "European Union" or "Euro Area" or "Euro zone")
	#         and map them to EU27. This catches labels like:
	#         "European Union - 27 Countries", "EURO AREA - 20 COUNTRIES", "Euro Zone (EA20)", etc.
	idx_eu_variants <- str_detect(x, "(?i)european\\s+union|euro\\s+area|euro\\s+zone|eurozone")
	out[idx_eu_variants] <- "EU27"

	# Step 2: Check for exact ISO3 codes.
	idx_iso3 <- str_detect(x, "^[A-Z]{3}$") & is.na(out)
	out[idx_iso3] <- x[idx_iso3]

	# Step 3: Check for exact EU matches (should be mostly covered by step 1, but be explicit).
	idx_eu <- x %in% c("EU", "EU27", "EU27_2020", "European Union") & is.na(out)
	out[idx_eu] <- "EU27"

	# Step 4: Check for ISO2 codes.
	idx_iso2 <- str_detect(x, "^[A-Z]{2}$") & is.na(out)
	out[idx_iso2] <- countrycode(x[idx_iso2], origin = "iso2c", destination = "iso3c")

	# Step 5: Check custom mappings.
	idx_other <- is.na(out)
	if (any(idx_other)) {
		out[idx_other] <- unname(custom[x[idx_other]])
	}

	# Step 6: Fall back to countrycode package for country names.
	idx_still_na <- is.na(out)
	if (any(idx_still_na)) {
		out[idx_still_na] <- countrycode(
			x[idx_still_na],
			origin = "country.name",
			destination = "iso3c",
			warn = FALSE
		)
	}

	out
}

# Convert monthly/annual labels to quarterly labels.
to_quarter <- function(x) {
	x <- str_trim(as.character(x))
	out <- rep(NA_character_, length(x))

	idx_q <- !is.na(x) & str_detect(x, "^\\d{4}-?Q[1-4]$")
	out[idx_q] <- str_remove(x[idx_q], "-")

	idx_q_nodash <- !is.na(x) & str_detect(x, "^\\d{4}Q[1-4]$")
	out[idx_q_nodash] <- x[idx_q_nodash]

	idx_m <- !is.na(x) & str_detect(x, "^\\d{4}-\\d{2}$")
	if (any(idx_m)) {
		dt <- suppressWarnings(ymd(paste0(x[idx_m], "-01")))
		out[idx_m] <- paste0(year(dt), "Q", quarter(dt))
	}

	idx_d <- !is.na(x) & str_detect(x, "^\\d{4}-\\d{2}-\\d{2}$")
	if (any(idx_d)) {
		dt <- suppressWarnings(ymd(x[idx_d]))
		out[idx_d] <- paste0(year(dt), "Q", quarter(dt))
	}

	idx_y <- !is.na(x) & str_detect(x, "^\\d{4}$")
	out[idx_y] <- paste0(x[idx_y], "Q4")

	out
}

# Aggregate multiple observations within a quarter.
aggregate_quarter <- function(x, method = "mean") {
	x <- as.numeric(x)
	x <- x[is.finite(x)]
	if (length(x) == 0) return(NA_real_)

	if (identical(method, "sum")) return(sum(x, na.rm = TRUE))
	if (identical(method, "last")) return(dplyr::last(x))
	if (identical(method, "first")) return(dplyr::first(x))
	mean(x, na.rm = TRUE)
}

# GDP-weighted mean that ignores missing values in either series.
weighted_mean_safe <- function(x, w) {
	ok <- is.finite(x) & is.finite(w) & !is.na(x) & !is.na(w)
	if (!any(ok)) return(NA_real_)
	sum(x[ok] * w[ok]) / sum(w[ok])
}

# EU27 member countries used to aggregate national OECD rates into an EU series.
eu27_iso3 <- c(
	"AUT", "BEL", "BGR", "HRV", "CYP", "CZE", "DNK", "EST", "FIN",
	"FRA", "DEU", "GRC", "HUN", "IRL", "ITA", "LVA", "LTU", "LUX",
	"MLT", "NLD", "POL", "PRT", "ROU", "SVK", "SVN", "ESP", "SWE"
)

ecb_policy_rate_q <- clean_ecb_policy_rate(ecb_rate)

# Standardize one raw panel into: country, date, indicator, value.
standardize_panel <- function(df, indicator, agg_method = "mean") {
	out <- df %>%
		mutate(
			iso3 = coalesce(to_iso3(country_raw), toupper(as.character(country_raw))),
			date = to_quarter(time_raw),
			value = as.numeric(value_raw)
		) %>%
		filter(!is.na(iso3), !is.na(date), is.finite(value)) %>%
		group_by(iso3, date) %>%
		summarise(value = aggregate_quarter(value, agg_method), .groups = "drop") %>%
		mutate(
			Pays = if_else(
				iso3 == "EU27",
				"European Union",
				coalesce(countrycode(iso3, origin = "iso3c", destination = "country.name", warn = FALSE), iso3)
			),
			indicator = indicator
		) %>%
		select(Pays, iso3, date, indicator, value)

	out %>% arrange(indicator, date, Pays)
}

# Generic OECD CSV reader with optional dimension filters.
read_oecd_long <- function(file, filters = list()) {
	df <- read_csv(file, show_col_types = FALSE)

	if (length(filters) > 0) {
		for (nm in names(filters)) {
			if (nm %in% names(df)) {
				df <- df %>% filter(.data[[nm]] %in% filters[[nm]])
			}
		}
	}

	df
}

# Clean one OECD long-format dataset.
clean_oecd_long <- function(file, indicator, filters = list(), agg_method = "mean") {
	df <- read_oecd_long(file, filters = filters)

	panel <- df %>%
		transmute(
			country_raw = REF_AREA,
			time_raw = TIME_PERIOD,
			value_raw = OBS_VALUE
		)

	standardize_panel(panel, indicator = indicator, agg_method = agg_method)
}

# Clean WDI deflator file and expand annual values to all quarters.
clean_wdi_deflator <- function(file, indicator = "gdp_deflator") {
	df <- read_csv(file, skip = 4, show_col_types = FALSE)

	year_cols <- names(df)[str_detect(names(df), "^\\d{4}$")]

	panel <- df %>%
		filter(`Indicator Code` == "NY.GDP.DEFL.ZS") %>%
		select(`Country Name`, `Country Code`, all_of(year_cols)) %>%
		pivot_longer(
			cols = all_of(year_cols),
			names_to = "time_raw",
			values_to = "value_raw"
		) %>%
		transmute(
			country_raw = `Country Code`,
			time_raw = time_raw,
			value_raw = value_raw
		)

	# Expand annual data to all 4 quarters
	panel <- panel %>%
		expand_grid(quarter = 1:4) %>%
		mutate(time_raw = paste0(time_raw, "Q", quarter)) %>%
		select(country_raw, time_raw, value_raw)

	standardize_panel(panel, indicator = indicator, agg_method = "mean")
}

# Clean Eurostat REER and keep one preferred REER definition per country/date.
clean_eurostat_reer <- function(file, indicator = "reer_ic37_41_cpi") {
	df <- read_csv(file, show_col_types = FALSE)

	# Normalize the source schema so we can handle small header variations.
	normalize_source_names <- function(df) {
		nms <- names(df)
		rename_map <- c(
			geo = "geo",
			TIME_PERIOD = "TIME_PERIOD",
			TIME_PERIOI = "TIME_PERIOD",
			TIME_PERIO = "TIME_PERIOD",
			exch_rt = "exch_rt",
			OBS_VALUE = "OBS_VALUE"
		)

		for (nm in names(rename_map)) {
			idx <- which(tolower(nms) == tolower(nm))
			if (length(idx) > 0) {
				nms[idx[1]] <- rename_map[[nm]]
			}
		}

		names(df) <- nms
		df
	}

	df <- normalize_source_names(df)

	if (!all(c("geo", "TIME_PERIOD", "OBS_VALUE") %in% names(df))) {
		stop(paste0(
			"eurostat_reer.csv must contain geo, TIME_PERIOD, and OBS_VALUE columns. Found: ",
			paste(names(df), collapse = ", ")
		))
	}

	if (!"exch_rt" %in% names(df)) {
		df$exch_rt <- NA_character_
	}

	panel <- df %>%
		mutate(exch_rt = as.character(exch_rt))

	# If the file contains Eurostat REER codes, keep the preferred definition.
	# Otherwise, fall back to the series present in the file so the panel is not empty.
	if ("exch_rt" %in% names(panel) && any(panel$exch_rt %in% c("REER_IC37_CPI", "REER_IC41_CPI", "REER_IC42_CPI"), na.rm = TRUE)) {
		panel <- panel %>%
			filter(exch_rt %in% c("REER_IC37_CPI", "REER_IC41_CPI", "REER_IC42_CPI")) %>%
			mutate(
				reer_priority = case_when(
					exch_rt == "REER_IC41_CPI" ~ 1L,
					exch_rt == "REER_IC42_CPI" ~ 2L,
					exch_rt == "REER_IC37_CPI" ~ 3L,
					TRUE ~ 99L
				)
			) %>%
			group_by(geo, TIME_PERIOD) %>%
			slice_min(order_by = reer_priority, n = 1, with_ties = FALSE) %>%
			ungroup()
	} else {
		panel <- panel %>%
			filter(!is.na(OBS_VALUE))
		message("Using exchange-rate series present in eurostat_reer.csv because REER codes were not found in exch_rt.")
	}

	panel <- panel %>%
		transmute(
			country_raw = geo,
			time_raw = TIME_PERIOD,
			value_raw = OBS_VALUE
		)

	standardize_panel(panel, indicator = indicator, agg_method = "mean")
}

# Dataset registry: one entry per indicator to clean.
datasets <- list(
	list(
		name = "oecd_real_gdp",
		indicator = "real_gdp_growth_qoq",
		cleaner = function() clean_oecd_long(
			file.path(raw_dir, "oecd_real_gdp.csv"),
			indicator = "real_gdp_growth_qoq",
			filters = list(
				TRANSACTION = "B1GQ",
				ADJUSTMENT = "Y",
				TRANSFORMATION = "G1"
			),
			agg_method = "first"
		)
	),
	list(
		name = "oecd_gdp_level",
		indicator = "gdp_level",
		cleaner = function() clean_oecd_long(
			file.path(raw_dir, "oecd_gdp_expenditures.csv"),
			indicator = "gdp_level",
			filters = list(
				TRANSACTION = "B1GQ",
				PRICE_BASE = "V",
				TRANSFORMATION = "N"
			),
			agg_method = "first"
		)
	),
	list(
		name = "oecd_gdp_expenditures",
		indicator = "final_consumption_expenditure",
		cleaner = function() clean_oecd_long(
			file.path(raw_dir, "oecd_gdp_expenditures.csv"),
			indicator = "final_consumption_expenditure",
			filters = list(
				TRANSACTION = "P3",
				PRICE_BASE = "V",
				TRANSFORMATION = c("N", "")
			),
			agg_method = "sum"
		)
	),
	list(
		name = "oecd_cpi",
		indicator = "cpi_index",
		cleaner = function() clean_oecd_long(
			file.path(raw_dir, "oecd_cpi.csv"),
			indicator = "cpi_index",
			filters = list(
				MEASURE = "CPI",
				EXPENDITURE = "_T",
				ADJUSTMENT = "N",
				TRANSFORMATION = "_Z"
			),
			agg_method = "last"
		)
	),
	list(
		name = "oecd_short_term_interest",
		indicator = "short_term_interest",
		cleaner = function() clean_oecd_long(
			file.path(raw_dir, "oecd_short_term_interest.csv"),
			indicator = "short_term_interest",
			filters = list(
				MEASURE = "IR3TIB",
				METHODOLOGY = "N"
			),
			agg_method = "mean"
		)
	),
	list(
		name = "oecd_trade",
		indicator = "trade_balance_usd_billion",
		cleaner = function() clean_oecd_long(
			file.path(raw_dir, "oecd_trade.csv"),
			indicator = "trade_balance_usd_billion",
			filters = list(
				TRADE_FLOW = "TB",
				COUNTERPART_AREA = "W",
				PRODUCT_TYPE = "C",
				ADJUSTMENT = "Y",
				TRANSFORMATION = "N"
			),
			agg_method = "sum"
		)
	),
	list(
		name = "oecd_deflator",
		indicator = "gdp_deflator",
		cleaner = function() clean_wdi_deflator(file.path(raw_dir, "oecd_deflator.csv"), "gdp_deflator")
	),
	list(
		name = "eurostat_reer",
		indicator = "reer_ic37_41_cpi",
		cleaner = function() clean_eurostat_reer(file.path(raw_dir, "eurostat_reer.csv"), "reer_ic37_41_cpi")
	)
)

# Run one dataset cleaning job with error isolation.
safe_run <- function(job) {
	message("Cleaning: ", job$name)
	res <- tryCatch(job$cleaner(), error = function(e) e)

	if (inherits(res, "error")) {
		message("  -> FAILED: ", res$message)
		return(NULL)
	}

	message("  -> OK: (", nrow(res), " rows)")
	res
}

# Execute all jobs and keep successful outputs only.
all_clean <- map(datasets, safe_run)

clean_ok <- all_clean[!map_lgl(all_clean, is.null)]

if (length(clean_ok) > 0) {
	# Diagnostic table: country/date coverage per indicator.
	coverage_report <- bind_rows(clean_ok) %>%
		group_by(indicator) %>%
		summarise(
			n_obs = n(),
			n_countries = n_distinct(iso3),
			countries = paste(sort(unique(iso3)), collapse = ", "),
			date_range = paste(min(date), "to", max(date)),
			.groups = "drop"
		)
	write_csv(coverage_report, file.path(out_dir, "coverage_report.csv"))
	message("\n=== COVERAGE REPORT ===")
	print(coverage_report)

	# Optional diagnostic output: sparse wide panel (can contain structural NAs).
	combined_sparse <- bind_rows(clean_ok) %>%
		select(Pays, iso3, date, indicator, value) %>%
		pivot_wider(
			names_from = indicator,
			values_from = value,
			values_fn = ~ dplyr::first(.x)
		)
	write_csv(combined_sparse, file.path(out_dir, "panel_all_indicators_sparse.csv"))
	message("Sparse panel (may contain NAs): ", file.path(out_dir, "panel_all_indicators_sparse.csv"), " (", nrow(combined_sparse), " rows)")

	# Final panel: keep the union of all available country-date observations.
	# This preserves source-provided aggregates (e.g., EU27) without manual reconstruction.
	eu_short_term_interest_q <- combined_sparse %>%
		filter(iso3 %in% eu27_iso3) %>%
		group_by(date) %>%
		summarise(
			eu_short_term_interest = weighted_mean_safe(short_term_interest, gdp_level),
			.groups = "drop"
		)

	panel_merged <- combined_sparse %>%
		left_join(ecb_policy_rate_q, by = "date") %>%
		left_join(eu_short_term_interest_q, by = "date")

	if ("short_term_interest" %in% names(panel_merged)) {
		panel_merged <- panel_merged %>%
			mutate(
				short_term_interest = if_else(
					iso3 == "EU27",
					coalesce(eu_short_term_interest, ecb_policy_rate),
					short_term_interest
				)
			)
	} else {
		panel_merged <- panel_merged %>%
			mutate(short_term_interest = if_else(iso3 == "EU27", coalesce(eu_short_term_interest, ecb_policy_rate), NA_real_))
	}

	panel_merged <- panel_merged %>%
		select(-ecb_policy_rate, -eu_short_term_interest) %>%
		arrange(date, Pays)

	write_csv(panel_merged, file.path(out_dir, "panel_all_indicators.csv"))
	message("Merged panel (with NAs where data missing): ", file.path(out_dir, "panel_all_indicators.csv"), " (", nrow(panel_merged), " rows, ", n_distinct(panel_merged$iso3), " countries)")
}

