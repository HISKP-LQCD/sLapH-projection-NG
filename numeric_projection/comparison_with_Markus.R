q_avail_actual <- names(all_prescriptions[[total_momentum_str]][[irrep]][['1']][['1']])

q_avail_actual_matrix <- apply(
  stringr::str_match(q_avail_actual, '(-?\\d)(-?\\d)(-?\\d)')[, 2:4],
  1:2,
  as.integer)

q_markus_avail_actual <- t(total_momentum_ref / 2 - t(q_avail_actual_matrix))

elements <- apply(
  q_markus_avail_actual,
  1,
  function (q) sprintf('p: %d, q: (%s), g: \\gamma_{5}, \\gamma_{5}', sum(total_momentum^2), paste0(sprintf('%.1f', q), collapse = ', ')))

df_avail_actual <- data.frame(str = q_avail_actual, element = elements, stringsAsFactors = FALSE)

operator_indices_path <- sprintf('reference/rho_p%d_%s_operator-indices.tsv', total_momentum_sq, irrep)
operator_indices <- read.table(operator_indices_path, sep = '\t', header = TRUE)

filtered <- operator_indices %>%
  filter(p_x == total_momentum[1],
         p_y == total_momentum[2],
         p_z == total_momentum[3],
         alpha == 1)
stopifnot(nrow(filtered) == 1)
operator_id <- filtered$id

gevp_indices_path <- sprintf('reference/rho_p%d_%s_gevp-indices.tsv', sum(total_momentum^2), irrep)
gevp_indices <- read.table(gevp_indices_path, sep = '\t', header = TRUE, stringsAsFactors = FALSE)

mapping <- left_join(df_avail_actual, gevp_indices, by = c('element')) %>%
  filter(!is.na(id))

correlator_matrix_indices <- expand.grid(q_source = mapping$str, q_sink = mapping$str, stringsAsFactors = FALSE)

correlator_matrix_indices %<>%
  tibble::as.tibble() %>%
  left_join(mapping, by = c('q_source' = 'str')) %>%
  rename(id_source = id) %>%
  left_join(mapping, by = c('q_sink' = 'str')) %>%
  rename(id_sink = id) %>%
  select(q_source, q_sink, id_source, id_sink)

correlator_matrix_indices %<>%
  mutate(markus_file_name = sprintf('reference/rho_p%d_%s_op%d_gevp%d.%d.tsv', sum(total_momentum^2), irrep, operator_id, id_source, id_sink))

load_target_config <- function (path) {
  target <- read.table(path, header = TRUE)
  
  target_single <- dplyr::filter(target, cnfg == config_number)
  target_single$value
}

correlator_matrix_indices$target <- mapply(load_target_config, correlator_matrix_indices$markus_file_name, SIMPLIFY = FALSE, USE.NAMES = FALSE)

load_actual_config <- function (q_source, q_sink) {
  irrep_col <- if (irrep %in% c('T1u', 'E')) '2' else '1'
  resolved[[total_momentum_str]][[irrep]][[irrep_col]][['1']][[q_source]][[q_sink]]
}

correlator_matrix_indices$actual <- mapply(load_actual_config, correlator_matrix_indices$q_source, correlator_matrix_indices$q_sink, SIMPLIFY = FALSE, USE.NAMES = FALSE)

correlator_matrix_indices$time <- list(1:length(correlator_matrix_indices$target[[1]]))

correlator_matrix <- correlator_matrix_indices %>%
  tidyr::unnest() %>%
  tidyr::gather(key, correlator, actual, target)

p <- ggplot(correlator_matrix, aes(x = time, y = abs(correlator), color = key)) +
  geom_point(position = position_dodge(width = 0.3)) +
  scale_y_log10() +
  facet_grid(q_source ~ q_sink) +
  labs(title = "Comparison between Martin's and Markus' projection code",
       subtitle = sprintf('P = %s, irrep = %s', total_momentum_str, irrep),
       x = 't',
       y = expression(abs(C(t))),
       color = 'Data')

dir.create('comparison')
ggsave(sprintf('comparison/comparison-%s-%s.pdf', total_momentum_str, irrep), p, width = 10, height = 10)
