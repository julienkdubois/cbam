suppressPackageStartupMessages({
	library(dplyr)
	library(readr)
	library(tidyr)
})

# Resolve script-relative paths so the script works from any working directory.
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

# Resolve the project data root by trying common locations.
resolve_data_root <- function() {
	candidates <- c(
		script_dir,
		file.path(script_dir, "data"),
		file.path(getwd(), "cbam", "data"),
		file.path(getwd(), "data")
	)

	for (p in candidates) {
		if (dir.exists(file.path(p, "clean"))) {
			return(normalizePath(p, winslash = "/", mustWork = FALSE))
		}
	}

	# Fallback requested by user: <cwd>/cbam/data
	normalizePath(file.path(getwd(), "cbam", "data"), winslash = "/", mustWork = FALSE)
}

data_root <- resolve_data_root()
data_dir <- file.path(data_root, "clean")
raw_dir <- file.path(data_root, "raw")
out_dir <- file.path(data_root, "clean")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Load the cleaned panel produced by cleaning.R.
panel <- read_csv(file.path(data_dir, "panel_all_indicators.csv"), show_col_types = FALSE)

message("Loaded panel: ", nrow(panel), " rows, ", n_distinct(panel$iso3), " countries")
message("Date range: ", min(panel$date), " to ", max(panel$date))
message("Variables: ", paste(setdiff(names(panel), c("Pays", "iso3", "date")), collapse = ", "))

# ===== 1) SELECT ESTIMATION WINDOW =====
# Edit these bounds to change the estimation sample.
start_date <- "2000Q1"
end_date <- "2025Q4"

panel_window <- panel %>%
	filter(date >= start_date, date <= end_date) %>%
	arrange(date, iso3)

# Sanity check: filtering by date should not drop columns.
panel_cols <- names(panel)
window_cols <- names(panel_window)
dropped_cols <- setdiff(panel_cols, window_cols)
if (length(dropped_cols) > 0) {
	message("Warning: columns dropped in panel_window: ", paste(dropped_cols, collapse = ", "))
}

# Normalize the expenditure column name in case readr repaired duplicates.
if (!"final_consumption_expenditure" %in% window_cols) {
	exp_candidates <- grep("^final_consumption_expenditure", window_cols, value = TRUE)
	if (length(exp_candidates) > 0) {
		panel_window <- panel_window %>% rename(final_consumption_expenditure = all_of(exp_candidates[1]))
		window_cols <- names(panel_window)
	}
}

if (!"final_consumption_expenditure" %in% names(panel_window)) {
	stop(paste0(
		"Missing final_consumption_expenditure in panel_window. Available columns: ",
		paste(names(panel_window), collapse = ", ")
	))
}

if (!"gdp_level" %in% names(panel_window)) {
	stop(paste0(
		"Missing gdp_level in panel_window. Available columns: ",
		paste(names(panel_window), collapse = ", ")
	))
}

message("\n=== ESTIMATION WINDOW ===")
message("Selected: ", start_date, " to ", end_date)
message("Rows in window: ", nrow(panel_window))
message("Unique dates: ", n_distinct(panel_window$date))

# Ensure critical columns exist even if they weren't produced by cleaning.R.
if (!"reer_ic37_41_cpi" %in% names(panel_window)) {
	message("Warning: reer_ic37_41_cpi column missing from panel. Creating with all NAs.")
	panel_window <- panel_window %>% mutate(reer_ic37_41_cpi = NA_real_)
}

if (!"real_gdp_growth_qoq" %in% names(panel_window)) {
	stop("Critical column real_gdp_growth_qoq is missing from panel.")
}

message("Available columns after NA-fill: ", paste(names(panel_window), collapse = ", "))

# ===== 2) VARIABLE TRANSFORMATIONS =====
# GDP growth is already in growth-rate form.
# Consumption is explicitly deflated by CPI before computing QoQ growth.
# CPI, REER, and other series are converted to QoQ percent changes.
# Note: NAs in source series are properly preserved through transformations.

# Helper: quarter-on-quarter percent change.
pct_change <- function(x) {
	((x - lag(x)) / lag(x)) * 100
}

# Helper: seasonal adjustment for quarterly series (fallback to raw values on failure).
seasonal_adjust_q <- function(x) {
	x_num <- as.numeric(x)
	n <- length(x_num)

	if (n < 8 || all(is.na(x_num))) {
		return(x_num)
	}

	idx <- seq_len(n)
	good <- which(!is.na(x_num))

	if (length(good) < 8) {
		return(x_num)
	}

	x_filled <- x_num
	if (length(good) < n) {
		x_filled[is.na(x_filled)] <- approx(
			x = good,
			y = x_num[good],
			xout = idx[is.na(x_filled)],
			rule = 2
		)$y
	}

	sa_fit <- tryCatch(
		stats::stl(stats::ts(x_filled, frequency = 4), s.window = "periodic", robust = TRUE),
		error = function(e) NULL
	)

	if (is.null(sa_fit)) {
		return(x_num)
	}

	sa <- as.numeric(x_filled - sa_fit$time.series[, "seasonal"])
	sa[is.na(x_num)] <- NA_real_
	sa
}

# Helper: convert a Date to a Dynare-friendly quarterly label.
to_quarter_label <- function(x) {
  x <- as.Date(x)
  month_num <- as.integer(format(x, "%m"))
	  paste0(format(x, "%Y"), "Q", ((month_num - 1L) %/% 3L) + 1L)
}

# FRED daily bilateral exchange rate -> quarterly average level.
read_fred_exchange_quarterly <- function(file) {
	df <- read_csv(file, show_col_types = FALSE)

	df %>%
		transmute(
			date = as.Date(observation_date),
			exchange_rate = as.numeric(DEXUSEU)
		) %>%
		filter(!is.na(date), is.finite(exchange_rate)) %>%
		mutate(date = to_quarter_label(date)) %>%
		group_by(date) %>%
		summarise(exchange_rate = mean(exchange_rate, na.rm = TRUE), .groups = "drop") %>%
		arrange(date)
}

# Census monthly bilateral trade flows -> quarterly trade balance in USD millions.
read_census_trade_quarterly <- function(file) {
	if (!requireNamespace("readxl", quietly = TRUE)) {
		stop("Package 'readxl' is required to read census_trade_balance.xlsx")
	}

	month_map <- c(
		IJAN = 1L, IFEB = 2L, IMAR = 3L, IAPR = 4L, IMAY = 5L, IJUN = 6L,
		IJUL = 7L, IAUG = 8L, ISEP = 9L, IOCT = 10L, INOV = 11L, IDEC = 12L,
		EJAN = 1L, EFEB = 2L, EMAR = 3L, EAPR = 4L, EMAY = 5L, EJUN = 6L,
		EJUL = 7L, EAUG = 8L, ESEP = 9L, EOCT = 10L, ENOV = 11L, EDEC = 12L
	)

	df <- readxl::read_excel(file, sheet = "country")
	month_cols <- intersect(names(df), names(month_map))

	if (length(month_cols) == 0) {
		stop("No monthly import/export columns found in census_trade_balance.xlsx")
	}

	df %>%
		pivot_longer(all_of(month_cols), names_to = "month_code", values_to = "value") %>%
		mutate(
			iso3 = case_when(
				grepl("European Union", CTYNAME, ignore.case = TRUE) ~ "EU27",
				grepl("United States", CTYNAME, ignore.case = TRUE) ~ "USA",
				TRUE ~ NA_character_
			),
			flow = substr(month_code, 1, 1),
			month_num = unname(month_map[month_code]),
			date = as.Date(sprintf("%04d-%02d-01", as.integer(year), month_num)),
			value = as.numeric(value)
		) %>%
		filter(!is.na(iso3), !is.na(date), is.finite(value)) %>%
		group_by(iso3, date, flow) %>%
		summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
		group_by(iso3, date = to_quarter_label(date)) %>%
		summarise(
			imports = sum(value[flow == "I"], na.rm = TRUE),
			exports = sum(value[flow == "E"], na.rm = TRUE),
			trade_balance_usd_million = exports - imports,
			.groups = "drop"
		) %>%
		arrange(date, iso3)
}

panel_transformed <- panel_window %>%
	group_by(iso3) %>%
	mutate(
		# Keep GDP growth as provided.
		gdp_growth = real_gdp_growth_qoq,
		
		# Seasonal adjustment for non-adjusted series.
		consumption_sa = seasonal_adjust_q(final_consumption_expenditure),
		cpi_sa = seasonal_adjust_q(cpi_index),

		# Explicitly deflate seasonally-adjusted consumption using seasonally-adjusted CPI.
		real_consumption = consumption_sa / cpi_sa,
		consumption_growth_real = pct_change(real_consumption),
		
		# Convert seasonally-adjusted CPI index to QoQ inflation.
		inflation_qoq = pct_change(cpi_sa),
		
		# Convert REER index to QoQ appreciation.
		reer_appreciation = pct_change(reer_ic37_41_cpi),
		
		# Keep short-term rate as provided.
		interest_rate = short_term_interest,
		
		# Keep trade balance in level and compute share of GDP.
		# trade_balance is in USD billions, while gdp_level is in millions.
		# Current convention: convert billions -> millions with *1000, then compute
		# TB as percent of GDP: 100 * (TB_millions / GDP_millions).
		trade_balance = trade_balance_usd_billion,
		trade_balance_share_gdp = 100 * (trade_balance_usd_billion * 1000) / gdp_level
	) %>%
	ungroup()

# Keep final columns used for estimation and aggregation.
dynare_data <- panel_transformed %>%
	select(
		date, iso3, Pays,
		gdp_growth,
		gdp_level,
		consumption_growth_real,
		inflation_qoq,
		interest_rate,
		reer_appreciation,
		trade_balance,
		trade_balance_share_gdp
	) %>%
	arrange(date, iso3)

# Drop first row per country because QoQ changes require a lag.
dynare_data <- dynare_data %>%
	group_by(iso3) %>%
	slice(-1) %>%
	ungroup()

message("\n=== TRANSFORMED DATA ===")
message("Rows after transformation: ", nrow(dynare_data))
message("Countries: ", n_distinct(dynare_data$iso3))
message("Date range: ", min(dynare_data$date), " to ", max(dynare_data$date))
message("\nVariables for Dynare:")
message("  - gdp_growth: real GDP growth QoQ (%)")
message("  - consumption_growth_real: SA consumption growth QoQ (%) after CPI deflation")
message("  - inflation_qoq: SA CPI inflation QoQ (%)")
message("  - interest_rate: short-term interest rate (%)")
message("  - reer_appreciation: REER appreciation QoQ (%)")
message("  - trade_balance: trade balance (USD billions)")
message("  - trade_balance_share_gdp: 100 * (trade_balance_usd_billion * 1000) / gdp_level")

# Print compact summary statistics for quick QA.
message("\n=== SUMMARY STATISTICS ===")
print(dynare_data %>%
	select(-date, -iso3, -Pays) %>%
	summarise(across(everything(), list(
		mean = ~ mean(., na.rm = TRUE),
		sd = ~ sd(., na.rm = TRUE),
		min = ~ min(., na.rm = TRUE),
		max = ~ max(., na.rm = TRUE),
		na_count = ~ sum(is.na(.))
	), .names = "{.col}_{.fn}"))
)

# Save final dataset for Dynare estimation.
out_file <- file.path(out_dir, "data_dynare.csv")
write_csv(dynare_data, out_file)
message("\n✓ Estimation data saved to: ", out_file)

# Optional visual QA checks (not required for export).
library(ggplot2)

# 1) Compare key series for major economies.
vars_to_plot <- c("gdp_growth","inflation_qoq","interest_rate","reer_appreciation")
focus_countries <- c("USA","GBR","JPN")

plot_data <- dynare_data %>%
  filter(iso3 %in% focus_countries) %>%
  pivot_longer(all_of(vars_to_plot))

ggplot(plot_data, aes(x=date, y=value, color=iso3, group=iso3)) +
  geom_line() +
  facet_wrap(~name, scales="free_y") +
  theme_minimal() +
	labs(title="Main macro series") %>%
	print()

# 2) Inspect tails and outliers by variable.
dynare_data %>%
  select(-date,-iso3,-Pays) %>%
  summarise(across(everything(), 
    list(p1=~quantile(.,0.01,na.rm=T), p99=~quantile(.,0.99,na.rm=T),
         min=~min(.,na.rm=T), max=~max(.,na.rm=T)))) %>%
  pivot_longer(everything()) %>% print(n=Inf)


# ===== 4) BILATERAL DYNARE PANEL (EU27 vs Rest of World) ===================
#
# Home    (h) : EU27
# Foreign (f) : RoW = Eurostat REER IC37/IC41 basket, ALL EU COUNTRIES EXCLUDED.
#
# Observable mapping (identical naming to bilateral EU27/USA dataset):
#   obs_dy_h / obs_dy_f  : real GDP growth QoQ (%)
#   obs_dc_h / obs_dc_f  : real consumption growth QoQ (%)
#   obs_pi_h / obs_pi_f  : CPI inflation QoQ (%)
#   obs_r_h  / obs_r_f   : short-term interest rate / 4 (quarterly)
#   obs_de               : EU27 REER index LEVEL (reer_ic37_41_cpi, base = 100)
#   obs_tb_h / obs_tb_f  : trade balance as % of own GDP
#
# Real exchange rate
#   The EU27 REER (IC37/IC41 CPI-based) is already a clean bilateral rate of
#   EU27 vs its non-EU trading partners — which is exactly the RoW definition
#   used here. It is used directly as a LEVEL index. No transformation, no
#   aggregation from the RoW side needed.
#
# Aggregation rules for RoW
#   GDP growth      : sum GDP levels across RoW countries, then QoQ % change
#   Consumption     : sum levels, deflate by GDP-weighted aggregate CPI, then QoQ %
#   Inflation       : QoQ % change of GDP-weighted aggregate CPI index
#   Interest rate   : GDP-weighted mean of country short-term rates
#   Trade balance   : sum of individual country TBs present in panel (% of RoW GDP)
# ---------------------------------------------------------------------------

# ---- 4.0  Define RoW country set (non-EU only) -----------------------------
eu_member_iso3 <- c(
  "AUT", "BEL", "BGR", "CYP", "CZE", "DEU", "DNK", "ESP", "EST",
  "FIN", "FRA", "GRC", "HRV", "HUN", "IRL", "ITA", "LTU", "LUX",
  "LVA", "MLT", "NLD", "POL", "PRT", "ROU", "SVK", "SVN", "SWE"
)

# Eurostat IC37/IC41 basket — non-EU entries only.
# SGP, THA, MYS are in the theoretical basket but absent from the panel.
reer_basket_all <- c(
  "USA", "JPN", "GBR", "CHE", "NOR", "ISL", "CAN", "AUS", "NZL",
  "CHN", "KOR", "HKG", "SGP", "IDN", "IND", "THA", "MYS",
  "BRA", "MEX", "ARG", "CHL", "COL", "CRI",
  "TUR", "RUS", "ZAF", "SAU", "ISR"
)

accidental_eu <- intersect(reer_basket_all, eu_member_iso3)
if (length(accidental_eu) > 0) {
  stop("RoW basket contains EU member(s): ", paste(accidental_eu, collapse = ", "))
}

row_countries_available <- intersect(reer_basket_all, unique(panel_window$iso3))

message("\n=== EU27 vs RoW SETUP ===")
message("RoW basket defined (non-EU) : ", length(reer_basket_all))
message("RoW countries in panel      : ", length(row_countries_available))
message("Missing from panel          : ",
        paste(setdiff(reer_basket_all, row_countries_available), collapse = ", "))

# Helper for plotting quarter labels as dates.
quarter_to_date <- function(q) {
	y <- as.integer(substr(q, 1, 4))
	qnum <- as.integer(sub(".*Q", "", q))
	month_num <- (qnum - 1L) * 3L + 1L
	as.Date(sprintf("%04d-%02d-01", y, month_num))
}

# EU27 home block used by the bilateral Dynare-style output.
eu_panel <- panel_transformed %>%
	filter(iso3 == "EU27") %>%
	arrange(date) %>%
	transmute(
		date,
		dy_h = gdp_growth,
		dc_h = consumption_growth_real,
		pi_h = inflation_qoq,
		r_h = interest_rate / 4,
		tb_full_h_bn = trade_balance,
		gdp_h = gdp_level
	)

if (nrow(eu_panel) == 0) {
	stop("EU27 block is missing from panel_transformed.")
}

de_level <- panel_transformed %>%
	filter(iso3 == "EU27") %>%
	arrange(date) %>%
	select(date, de = reer_ic37_41_cpi)

# ---- 4.1  RoW aggregate — level aggregation first (robust method) ----------
# On agrège d'abord les niveaux au niveau RoW, puis on calcule les croissances.
# Cela évite les pics artificiels liés à la moyenne de taux de croissance pays par pays.

row_panel <- panel_transformed %>%
  filter(iso3 %in% row_countries_available)

row_levels <- row_panel %>%
  group_by(date) %>%
  summarise(
		# Niveaux agrégés pour les variables qui supportent une somme.
		gdp_level_row     = sum(gdp_level, na.rm = TRUE),
		consumption_nominal_row = sum(consumption_sa, na.rm = TRUE),
		cpi_row           = weighted.mean(cpi_sa, gdp_level, na.rm = TRUE),
    interest_row_wavg  = weighted.mean(short_term_interest, gdp_level, na.rm = TRUE),

		# Sommes pour les ratios de balance commerciale.
    tb_row_usd_bn     = sum(trade_balance_usd_billion, na.rm = TRUE),
    
    n_countries       = sum(!is.na(gdp_level)),
    gdp_coverage_frac = sum(gdp_level[!is.na(gdp_growth)], na.rm = TRUE) / sum(gdp_level, na.rm = TRUE),
		consumption_coverage_frac = sum(gdp_level[!is.na(consumption_sa) & !is.na(cpi_sa)], na.rm = TRUE) / sum(gdp_level, na.rm = TRUE),
    .groups = "drop"
  ) %>%
	arrange(date) %>%
	mutate(
		dy_f = pct_change(gdp_level_row),
		dc_f = pct_change(consumption_nominal_row / cpi_row),
		pi_f = pct_change(cpi_row),
		r_f = interest_row_wavg / 4,
		tb_f_usd_bn = tb_row_usd_bn,
		gdp_f = gdp_level_row
	)

# Guardrail: banish isolated spikes in RoW consumption growth and interpolate them.
# This targets discontinuities like a one-quarter +120% jump caused by bad source data.
spike_threshold <- 120
dc_f_spike_dates <- row_levels %>%
	filter(is.finite(dc_f) & abs(dc_f) > spike_threshold) %>%
	pull(date)

if (length(dc_f_spike_dates) > 0) {
	message(
		"Warning: removing dc_f spike(s) above ", spike_threshold, "% on: ",
		paste(dc_f_spike_dates, collapse = ", ")
	)

	row_levels <- row_levels %>%
		mutate(dc_f = if_else(is.finite(dc_f) & abs(dc_f) > spike_threshold, NA_real_, dc_f)) %>%
		mutate(
			dc_f = {
				valid_idx <- which(!is.na(dc_f))
				if (length(valid_idx) < 2) {
					dc_f
				} else {
					stats::approx(
						x = valid_idx,
						y = dc_f[valid_idx],
						xout = seq_along(dc_f),
						rule = 2
					)$y
				}
			}
		)
}

## Robust outlier cleaning across key RoW series
# Helper: replace extreme outliers (median +/- mult * MAD) with linear interpolation
clean_outliers_mad <- function(x, mult = 5) {
	if (sum(!is.na(x)) < 3) return(x)
	med <- median(x, na.rm = TRUE)
	m <- mad(x, constant = 1, na.rm = TRUE)
	if (!is.finite(m) || m == 0) return(x)
	out_idx <- which(!is.na(x) & abs(x - med) > mult * m)
	if (length(out_idx) == 0) return(x)
	x[out_idx] <- NA_real_
	valid_idx <- which(!is.na(x))
	if (length(valid_idx) < 2) return(x)
	stats::approx(x = valid_idx, y = x[valid_idx], xout = seq_along(x), rule = 2)$y
}

# Apply cleaning to several series and log removed dates
for (series_name in c("dy_f", "dc_f", "pi_f", "r_f")) {
	if (!series_name %in% names(row_levels)) next
	series_vec <- row_levels[[series_name]]
	med <- median(series_vec, na.rm = TRUE)
	m <- mad(series_vec, constant = 1, na.rm = TRUE)
	if (!is.finite(m) || m == 0) next
	out_idx <- which(!is.na(series_vec) & abs(series_vec - med) > 5 * m)
	if (length(out_idx) > 0) {
		msg_dates <- paste(row_levels$date[out_idx], collapse = ", ")
		message(sprintf("Info: cleaning %s outliers on: %s", series_name, msg_dates))
		row_levels[[series_name]] <- clean_outliers_mad(series_vec, mult = 5)
	}
}

low_consumption_coverage_dates <- row_levels %>%
	filter(consumption_coverage_frac < 0.9) %>%
 	pull(date)

if (length(low_consumption_coverage_dates) > 0) {
	message(
		"Warning: RoW consumption coverage below 90% on: ",
		paste(low_consumption_coverage_dates, collapse = ", ")
	)
}

# ---- 4.2  RoW foreign block formatting -------------------------------------
row_block <- row_levels %>%
  transmute(
    date,
		dy_f,
		dc_f,
		pi_f,
		r_f,
		tb_f_usd_bn,
		gdp_f,
    n_countries,
    gdp_coverage_frac
  )

# ---- 4.6  Join all blocks and final calculations ---------------------------
dynare_eu_row <- eu_panel %>%
  inner_join(row_block,  by = "date") %>%
  inner_join(de_level,   by = "date") %>%
  mutate(
    # Calcul des ratios de balance commerciale
    tb_h = 100 * (tb_full_h_bn * 1000) / gdp_h,
    tb_f = 100 * (tb_f_usd_bn  * 1000) / gdp_f,
    # Calcul indispensable : variation du REER (car Dynare attend une croissance)
    obs_de = 100 * (de / lag(de) - 1)
  ) %>%
  arrange(date)

# ---- 4.7  Nettoyage robuste et Troncature automatique -----------------------

# 1. Liste des variables critiques pour Dynare
core_obs_eu_row <- c("dy_h", "dy_f", "dc_h", "dc_f", "pi_h", "pi_f", "r_h", "r_f", "obs_de", "tb_h", "tb_f")

# 2. Interpolation des trous au milieu des séries (middle NAs)
# Si un pays manque une donnée un trimestre, weighted.mean peut renvoyer NA.
# On remplit ces trous par interpolation linéaire pour ne pas perdre toute la série.
dynare_eu_row <- dynare_eu_row %>%
  mutate(across(all_of(core_obs_eu_row), function(x) {
    if (sum(!is.na(x)) < 2) return(x) # Trop peu de données
    valid_idx <- which(!is.na(x))
    # On remplit uniquement les NA situés ENTRE des valeurs valides
    approx(x = valid_idx, y = x[valid_idx], xout = seq_along(x), rule = 1)$y
  }))

# 3. Identification de la ligne de début et de fin (sans aucun NA)
dynare_eu_row$complete_row <- complete.cases(dynare_eu_row[, core_obs_eu_row])

if (!any(dynare_eu_row$complete_row)) {
  # DEBUG : Si aucune ligne n'est complète, on affiche quelle colonne est vide
  na_counts <- colSums(is.na(dynare_eu_row[, core_obs_eu_row]))
  print(na_counts)
  stop("Aucune ligne n'est complète. Vérifiez les colonnes ci-dessus qui ont trop de NA.")
}

first_idx <- min(which(dynare_eu_row$complete_row))
last_idx  <- max(which(dynare_eu_row$complete_row))

# 4. Troncature finale
dynare_eu_row_final <- dynare_eu_row[first_idx:last_idx, ] %>%
  select(-complete_row)

# 5. Vérification ultime et diagnostic
missing_after <- colSums(is.na(dynare_eu_row_final[, core_obs_eu_row]))
if (any(missing_after > 0)) {
  message("Erreur critique : Des NA persistent dans ces colonnes :")
  print(missing_after[missing_after > 0])
  stop("L'interpolation a échoué à boucher les trous internes.")
}

message("✓ Données prêtes. Début : ", dynare_eu_row_final$date[1], " Fin : ", tail(dynare_eu_row_final$date, 1))
# ---- 4.8  Dynare-ready observables (obs_ naming convention) ----------------
dynare_eu_row_obs <- dynare_eu_row_final %>%
  transmute(
    date,
    obs_dy_h = dy_h,
    obs_dy_f = dy_f,
    obs_dc_h = dc_h,
    obs_dc_f = dc_f,
    obs_pi_h = pi_h,
    obs_pi_f = pi_f,
    obs_r_h  = r_h,
    obs_r_f  = r_f,
		obs_de = obs_de,
    obs_tb_h = tb_h,
    obs_tb_f = tb_f
  )

# ---- 4.9  Save outputs -----------------------------------------------------
# Save the Dynare-ready data with the date column first so dseries can parse it.
dynare_data_final <- dynare_eu_row_obs

write_csv(dynare_data_final,
		  file.path(out_dir, "dynare_data_eu_row.csv"),
		  col_names = TRUE)

# ---- 4.10  QA plots --------------------------------------------------------
macro_plot_data <- dynare_eu_row_final %>%
  transmute(
    date = quarter_to_date(date),
    dy_h, dy_f, dc_h, dc_f,
    pi_h, pi_f, r_h, r_f,
    tb_h, tb_f
  ) %>%
  pivot_longer(
    cols          = -date,
    names_to      = c("variable", "zone"),
    names_pattern = "(.+)_(h|f)$",
    values_to     = "value"
  ) %>%
  mutate(zone = dplyr::recode(zone, h = "EU27", f = "RoW"))

p_macro <- ggplot(macro_plot_data,
                  aes(x = date, y = value, color = zone, group = zone)) +
  geom_line(linewidth = 0.5) +
  facet_wrap(~variable, scales = "free_y") +
  theme_minimal() +
  labs(title = "EU27 (home) vs Rest of World (foreign) — macro observables",
       x = "Date", y = NULL, color = "Area")

print(p_macro)
ggsave(file.path(out_dir, "dynare_bilateral_eu_row_macro_plot.png"),
       p_macro, width = 12, height = 7, dpi = 150)

p_reer <- dynare_eu_row_final %>%
  transmute(date = quarter_to_date(date), `EU27 REER (obs_de)` = de) %>%
  pivot_longer(-date) %>%
  ggplot(aes(x = date, y = value, color = name)) +
  geom_line(linewidth = 0.6) +
  theme_minimal() +
  labs(title = "Real exchange rate — EU27 REER index level (obs_de)",
       subtitle = "Source: Eurostat IC37/IC41 CPI-based REER. Base = 100 (reference year).",
       x = "Date", y = "Index", color = NULL)

print(p_reer)
ggsave(file.path(out_dir, "dynare_bilateral_eu_row_reer_plot.png"),
       p_reer, width = 10, height = 5, dpi = 150)

message("✓ QA plots saved to : ", out_dir)

# # ===== 3) BILATERAL DYNARE PANEL (EU27 vs USA) =====
# # Build final observables for a two-country setup.
# dynare_start_date <- "2002Q1"
# 
# fred_exchange_q <- read_fred_exchange_quarterly(file.path(raw_dir, "fred_exchange_rate.csv"))
# 
# cpi_bilateral_q <- panel_transformed %>%
# 	filter(iso3 %in% c("EU27", "USA")) %>%
# 	select(date, iso3, cpi_sa) %>%
# 	pivot_wider(names_from = iso3, values_from = cpi_sa) %>%
# 	rename(cpi_h = EU27, cpi_f = USA)
# 
# census_tb_q <- read_census_trade_quarterly(file.path(raw_dir, "census_trade_balance.xlsx")) %>%
# 	filter(iso3 == "EU27") %>%
# 	left_join(
# 		panel_window %>%
# 			filter(iso3 == "USA") %>%
# 			select(date, gdp_level_usa = gdp_level),
# 		by = "date"
# 	) %>%
# 	left_join(
# 		panel_window %>%
# 			filter(iso3 == "EU27") %>%
# 			select(date, gdp_level_eu = gdp_level),
# 		by = "date"
# 	) %>%
# 	mutate(
# 		# Census EU27 row gives bilateral US balance vis-a-vis EU: exports - imports.
# 		# We define foreign (US) TB directly, and home (EU) as the opposite sign.
# 		tb_f = 100 * trade_balance_usd_million / gdp_level_usa,
# 		tb_h = -100 * trade_balance_usd_million / gdp_level_eu
# 	) %>%
# 	select(date, tb_h, tb_f)
# 
# home_block <- dynare_data %>%
# 	filter(iso3 == "EU27") %>%
# 	transmute(
# 		date,
# 		dy_h = gdp_growth,
# 		dc_h = consumption_growth_real,
# 		pi_h = inflation_qoq,
# 		r_h = interest_rate / 4,
# 		reer_h = reer_appreciation
# 	)
# 
# foreign_block <- dynare_data %>%
# 	filter(iso3 == "USA") %>%
# 	transmute(
# 		date,
# 		dy_f = gdp_growth,
# 		dc_f = consumption_growth_real,
# 		pi_f = inflation_qoq,
# 		r_f = interest_rate / 4,
# 		reer_f = reer_appreciation
# 	)
# 
# dynare_bilateral <- home_block %>%
# 	inner_join(foreign_block, by = "date") %>%
# 	left_join(fred_exchange_q, by = "date") %>%
# 	left_join(cpi_bilateral_q, by = "date") %>%
# 	left_join(census_tb_q, by = "date") %>%
# 	mutate(
# 		real_exchange_rate = exchange_rate * (cpi_f / cpi_h),
# 			de = real_exchange_rate
# 	) %>%
# 	select(date, dy_h, dy_f, dc_h, dc_f, pi_h, pi_f, r_h, r_f, de, tb_h, tb_f) %>%
# 	arrange(date)
# 
# # Keep the latest possible end date such that the full sample has no NA.
# dynare_bilateral <- dynare_bilateral %>%
# 	filter(date >= dynare_start_date)
# 
# if (nrow(dynare_bilateral) == 0) {
# 	stop("No bilateral observations available from 2002Q1.")
# }
# 
# dynare_bilateral <- dynare_bilateral %>%
# 	mutate(complete_row = if_all(-date, ~ !is.na(.)))
# 
# # Find first and last complete observations
# first_complete_idx <- which(dynare_bilateral$complete_row)[1]
# last_complete_idx <- tail(which(dynare_bilateral$complete_row), 1)
# 
# if (is.na(first_complete_idx) || is.na(last_complete_idx)) {
# 	stop("No complete (non-NA) bilateral observations found in 2002Q1 onwards.")
# }
# 
# if (first_complete_idx > 1) {
# 	first_incomplete <- dynare_bilateral$date[1:(first_complete_idx - 1)]
# 	message("Warning: Bilateral panel has NAs at start (", paste(first_incomplete, collapse = ", "), ")")
# 	message("         Starting from first complete observation: ", dynare_bilateral$date[first_complete_idx])
# }
# 
# dynare_start_date_actual <- dynare_bilateral$date[first_complete_idx]
# dynare_end_date <- dynare_bilateral$date[last_complete_idx]
# 
# dynare_data <- dynare_bilateral %>%
# 	filter(date <= dynare_end_date) %>%
# 	select(-complete_row)
# 
# if (any(!complete.cases(dynare_data %>% select(-date)))) {
# 	stop("NA values remain in final dynare_data despite end-date truncation.")
# }
# 
# message("\n=== BILATERAL DYNARE DATA ===")
# message("Countries: EU27 (home) vs USA (foreign)")
# message("Period: ", dynare_start_date_actual, " to ", dynare_end_date)
# message("Rows: ", nrow(dynare_data))
# 
# raw_dynare_data <- dynare_data
# 
# dynare_observables <- raw_dynare_data %>%
# 	rename(
# 		obs_dy_h = dy_h,
# 		obs_dy_f = dy_f,
# 		obs_dc_h = dc_h,
# 		obs_dc_f = dc_f,
# 		obs_pi_h = pi_h,
# 		obs_pi_f = pi_f,
# 		obs_r_h  = r_h,
# 		obs_r_f  = r_f,
# 		obs_de   = de,
# 		obs_tb_h = tb_h,
# 		obs_tb_f = tb_f
# 	)
# 
# dynare_data <- dynare_observables %>%
# 	select(-date)
# 
# legacy_out_file <- file.path(out_dir, "data_dynare.csv")
# write_csv(raw_dynare_data, legacy_out_file)
# message("✓ Bilateral raw data saved to: ", legacy_out_file)
# 
# dynare_out_file <- file.path(out_dir, "dynare_data_bilateral_2002Q1_2019Q4.csv")
# write_csv(dynare_data, dynare_out_file)
# message("✓ Dynare-ready data saved to: ", dynare_out_file)
# # 
# # 3) Plot saved bilateral EU27/USA series (the same columns saved in data_dynare.csv).
# quarter_to_date <- function(q) {
# 	y <- substr(q, 1, 4)
# 	qnum <- as.integer(sub(".*Q", "", q))
# 	m <- (qnum - 1L) * 3L + 1L
# 	as.Date(sprintf("%s-%02d-01", y, m))
# }
# 
# bilat_plot_pairs <- raw_dynare_data %>%
# 	transmute(
# 		date = quarter_to_date(date),
# 		dy_h, dy_f,
# 		dc_h, dc_f,
# 		pi_h, pi_f,
# 		r_h, r_f,
# 		tb_h, tb_f
# 	) %>%
# 	pivot_longer(
# 		cols = -date,
# 		names_to = c("variable", "zone"),
# 		names_pattern = "(.*)_(h|f)",
# 		values_to = "value"
# 	) %>%
# 	mutate(zone = dplyr::recode(zone, h = "EU27", f = "USA"))
# 
# de_plot_data <- raw_dynare_data %>%
# 	transmute(
# 		date = quarter_to_date(date),
# 		variable = "de",
# 		zone = "Bilateral",
# 		value = de
# 	)
# 
# bilat_plot_data <- bind_rows(bilat_plot_pairs, de_plot_data)
# 
# bilat_plot <- ggplot(bilat_plot_data, aes(x = date, y = value, color = zone, group = zone)) +
# 	geom_line(linewidth = 0.5) +
# 	facet_wrap(~variable, scales = "free_y") +
# 	theme_minimal() +
# 	labs(title = "Series saved for EU27 and USA", x = "Date", y = "Value", color = "Area")
# 
# print(bilat_plot)
# 
# plot_file <- file.path(out_dir, "dynare_bilateral_us_eu_plot.png")
# ggsave(filename = plot_file, plot = bilat_plot, width = 12, height = 7, dpi = 150)
# message("✓ Bilateral plot saved to: ", plot_file)
