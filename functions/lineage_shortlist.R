# generalize lineages for color schemes in lineage frequency plot

lineage_shortlist = function(df) {
  
  df = df %>%
    mutate(ShortLin = pango_lineage)
  
  df_filt = df %>%
  mutate(ShortLin = case_when(
    endsWith(ShortLin, "2.12.1") ~ "BA.2.12.1",
    startsWith(ShortLin, "BG") ~ "BA.2.12.1",
    startsWith(ShortLin, "BA.2.75") ~ "BA.2.75",
    startsWith(ShortLin, "BA.1") ~ "BA.1",
    startsWith(ShortLin, "BA.2") ~ "BA.2",
    startsWith(ShortLin, "BK") ~ "BA.2",
    startsWith(ShortLin, "BP") ~ "BA.2",
    startsWith(ShortLin, "BL") ~ "BA.2.75",
    startsWith(ShortLin, "BH") ~ "BA.2",
    startsWith(ShortLin, "BM") ~ "BA.2",
    startsWith(ShortLin, "CH") ~ "BA.2",
    startsWith(ShortLin, "BJ") ~ "BA.2",
    startsWith(ShortLin, "BR") ~ "BA.2",
    startsWith(ShortLin, "DS") ~ "BA.2",
    startsWith(ShortLin, "BY") ~ "BA.2",
    startsWith(ShortLin, "DV") ~ "BA.2.75",
    startsWith(ShortLin, "BA.4") ~ "BA.4",
    startsWith(ShortLin, "BA.4.6") ~ "BA.4.6",
    startsWith(ShortLin, "BA.5") ~ "BA.5",
    startsWith(ShortLin, "BU") ~ "BA.5",
    startsWith(ShortLin, "BV") ~ "BA.5",
    startsWith(ShortLin, "BE") ~ "BA.5",
    startsWith(ShortLin, "BT") ~ "BA.5",
    startsWith(ShortLin, "BF") ~ "BA.5",
    startsWith(ShortLin, "BQ.1") ~ "BA.5",
    startsWith(ShortLin, "BQ.2") ~ "BA.5",
    startsWith(ShortLin, "BQ.3") ~ "BA.5",
    startsWith(ShortLin, "CN") ~ "BA.5",
    startsWith(ShortLin, "CD") ~ "BA.5",
    startsWith(ShortLin, "CE") ~ "BA.5",
    startsWith(ShortLin, "CF") ~ "BA.5",
    startsWith(ShortLin, "CL") ~ "BA.5",
    startsWith(ShortLin, "CG") ~ "BA.5",
    startsWith(ShortLin, "CK") ~ "BA.5",
    startsWith(ShortLin, "BZ") ~ "BA.5",
    startsWith(ShortLin, "XBB") ~ "XBB",
    TRUE ~ "Other"
  ))
  
  return(df_filt)

}