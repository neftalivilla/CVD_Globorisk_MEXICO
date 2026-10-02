## Data Analysis: Neftali Eduardo Antonio Villa & Omar Yaxmehen Bello-Chavolla
## Latest version of Analysis: April 2026
## For any question regarding analysis contact:
## Neftali Eduardo Antonio Villa: eduardo.villa@cardiologia.org.mx

#####Dataset: Library#####

library(dplyr)
library(tidyr)
library(stringr)
library(globorisk)
library(purrr)
library(flextable)
library(officer)

#####Dataset: Working Directory#####

setwd('/Users/nefoantonio/Library/CloudStorage/GoogleDrive-nantonio@facmed.unam.mx/.shortcut-targets-by-id/1RG-1Sg0NBmkz3AlzY5XL6XJ0mezx6oNg/Mexico City Cohort Study/Proyectos/0_Finalized projects/CVD risk scores/Table_Risks')

#####Dataset: Functions#####

globorisk_recalib_fast <- function(sex, age, sbp, tcol_mmol, diab, smoke, validate = TRUE) {
  
  if (validate) {
    # Verificación mínima y vectorizada
    vars <- list(sex = sex, age = age, sbp = sbp,
                 tcol_mmol = tcol_mmol, diab = diab, smoke = smoke)
    
    no_numericas <- names(vars)[!vapply(vars, is.numeric, logical(1))]
    if (length(no_numericas) > 0) {
      stop(sprintf(
        "Estas variables deben ser numéricas: %s",
        paste(no_numericas, collapse = ", ")
      ))
    }
    
    if (anyNA(sex) || anyNA(age) || anyNA(sbp) || anyNA(tcol_mmol) || anyNA(diab) || anyNA(smoke)) {
      stop("No se permiten valores NA en las variables de entrada.")
    }
    
    if (any(sex != 0 & sex != 1)) {
      stop("sexo debe ser 0 = hombre y 1 = mujer")
    }
    
    if (any(diab != 0 & diab != 1)) {
      stop("diabetes debe ser 0 = no y 1 = sí")
    }
    
    if (any(smoke != 0 & smoke != 1)) {
      stop("tabaquismo debe ser 0 = no y 1 = sí")
    }
    
    n <- length(sex)
    if (!all(lengths(vars) == n)) {
      stop("Todas las variables deben tener la misma longitud.")
    }
  }
  
  # Riesgo base Globorisk
  score <- globorisk(
    sex = sex,
    age = age,
    sbp = sbp,
    tc = tcol_mmol,
    dm = diab,
    smk = smoke,
    iso = "MEX",
    year = 2020,
    version = "fatal",
    type = "risk"
  )
  
  # Parámetros por sexo: 0 = hombre, 1 = mujer
  alpha <- c(-0.8639137, -0.881977)
  beta  <- c( 0.7745622,  1.002298)
  
  surv_model <- 1 - score
  log_term <- log(-log(surv_model))
  
  # sex + 1 sirve para indexar: 0->1 (hombre), 1->2 (mujer)
  log_recalib <- alpha[sex + 1] + beta[sex + 1] * log_term
  surv_recalib <- exp(-exp(log_recalib))
  rsk_recalib <- 1 - surv_recalib
  
  round(rsk_recalib * 100, 0)
}
make_globorisk_matrix <- function(data, sex_value, diab_value, smoke_value) {
  
  data %>%
    filter(
      sex == sex_value,
      diab == diab_value,
      smoke == smoke_value
    ) %>%
    select(age_group, sbp, tcol_mmol, risk_recalib) %>%
    mutate(risk_recalib = round(risk_recalib, 1)) %>%
    pivot_wider(
      names_from = tcol_mmol,
      values_from = risk_recalib
    ) %>%
    arrange(age_group, sbp) %>%
    rename(
      `Edad` = age_group,
      `PAS` = sbp
    )
}
make_ft <- function(tabla, titulo = NULL) {
  
  ft <- flextable(tabla)
  
  ft <- merge_v(ft, j = "Edad")
  ft <- valign(ft, j = "Edad", valign = "center")
  ft <- align(ft, align = "center", part = "all")
  ft <- autofit(ft)
  ft <- theme_booktabs(ft)
  
  if (!is.null(titulo)) {
    ft <- set_caption(ft, caption = titulo)
  }
  
  ft
}
color_globorisk <- function(ft, data_cols) {
  
  # cortes similares a la figura
  ft <- bg(ft, j = data_cols, i = ~ get(data_cols[1]) < 5, bg = "#00B050")
  ft <- bg(ft, j = data_cols, i = ~ get(data_cols[1]) >= 5  & get(data_cols[1]) < 10, bg = "#92D050")
  ft <- bg(ft, j = data_cols, i = ~ get(data_cols[1]) >= 10 & get(data_cols[1]) < 20, bg = "#FFFF00")
  ft <- bg(ft, j = data_cols, i = ~ get(data_cols[1]) >= 20 & get(data_cols[1]) < 30, bg = "#FFC000")
  ft <- bg(ft, j = data_cols, i = ~ get(data_cols[1]) >= 30 & get(data_cols[1]) < 40, bg = "#FF9999")
  ft <- bg(ft, j = data_cols, i = ~ get(data_cols[1]) >= 40, bg = "#FF0000")
  
  ft
}

apply_risk_colors <- function(ft, tabla) {
  
  risk_cols <- c("3", "4", "5", "6", "7")
  
  # Paleta extendida: cortes de 5 en 5 a 8 niveles
  pal <- list(
    p1 = "#4CAF50",  # < 5
    p2 = "#8BC34A",  # 5–9
    p3 = "#CDDC39",  # 10–14
    p4 = "#FFEB3B",  # 15–19
    p5 = "#FFC107",  # 20–24
    p6 = "#FF9800",  # 25–29
    p7 = "#F44336",  # 30–39
    p8 = "#7B0000"   # ≥ 40
  )
  
  for (col in risk_cols) {
    vals <- suppressWarnings(as.numeric(tabla[[col]]))
    vals <- pmin(vals, 40)
    
    ft <- bg(ft, j = col, i = which(!is.na(vals) & vals < 5),                   bg = pal$p1)
    ft <- bg(ft, j = col, i = which(!is.na(vals) & vals >= 5  & vals < 10),     bg = pal$p2)
    ft <- bg(ft, j = col, i = which(!is.na(vals) & vals >= 10 & vals < 15),     bg = pal$p3)
    ft <- bg(ft, j = col, i = which(!is.na(vals) & vals >= 15 & vals < 20),     bg = pal$p4)
    ft <- bg(ft, j = col, i = which(!is.na(vals) & vals >= 20 & vals < 25),     bg = pal$p5)
    ft <- bg(ft, j = col, i = which(!is.na(vals) & vals >= 25 & vals < 30),     bg = pal$p6)
    ft <- bg(ft, j = col, i = which(!is.na(vals) & vals >= 30 & vals < 40),     bg = pal$p7)
    ft <- bg(ft, j = col, i = which(!is.na(vals) & vals >= 40),                 bg = pal$p8)
  }
  
  ft
}

#####Dataset: Dummy Variables for Risk Prediction (in mmol/dl)#####

# 1) Definir los niveles de cada dimensión de la tabla
sex <- c(0, 1)              # 0 = hombre, 1 = mujer
diab <- c(0, 1)             # 0 = no diabético, 1 = diabético
smoke <- c(0, 1)            # 0 = no fumador, 1 = fumador
sbp <- c(120, 140, 160, 180)
tcol_mmol <- c(3, 4, 5, 6, 7)

# 2) Definir el mapeo de grupos de edad a una edad numérica representativa
age_map <- tibble(
  age_group = c("40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74"),
  age = c(44, 49, 54, 59, 64, 69, 74))

# 3) Generar todas las permutaciones posibles de las celdas
globorisk_dummy <- tidyr::crossing(
  sex = sex,
  diab = diab,
  smoke = smoke,
  age_map,
  sbp = sbp,
  tcol_mmol = tcol_mmol) %>%
  mutate(cell_id = row_number()) %>%
  arrange(
    desc(sex),          # Mujeres primero, luego Hombres
    diab,               # No diabéticos primero
    smoke,              # No fumadores primero
    desc(age),          # 70-74 arriba
    desc(sbp),          # 180 arriba
    tcol_mmol           # 3 a 7 de izquierda a derecha
  ) %>%
  select(
    cell_id,
    sex,
    age_group, age,
    sbp,
    tcol_mmol,
    diab, 
    smoke)


#####Dataset: Dummy Variables for Risk Prediction (in mg/dl)#####

# 1) Definir los niveles de cada dimensión de la tabla
sex <- c(0, 1)              # 0 = hombre, 1 = mujer
diab <- c(0, 1)             # 0 = no diabético, 1 = diabético
smoke <- c(0, 1)            # 0 = no fumador, 1 = fumador
sbp <- c(120, 140, 160, 180)
tcol_mmol <- c(2.6, 3.9, 5.2, 6.5, 7.8)

# 2) Definir el mapeo de grupos de edad a una edad numérica representativa
age_map <- tibble(
  age_group = c("40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74"),
  age = c(44, 49, 54, 59, 64, 69, 74))

# 3) Generar todas las permutaciones posibles de las celdas
globorisk_dummy_mg_dl <- tidyr::crossing(
  sex = sex,
  diab = diab,
  smoke = smoke,
  age_map,
  sbp = sbp,
  tcol_mmol = tcol_mmol) %>%
  mutate(cell_id = row_number()) %>%
  arrange(
    desc(sex),          # Mujeres primero, luego Hombres
    diab,               # No diabéticos primero
    smoke,              # No fumadores primero
    desc(age),          # 70-74 arriba
    desc(sbp),          # 180 arriba
    tcol_mmol           # 3 a 7 de izquierda a derecha
  ) %>%
  select(
    cell_id,
    sex,
    age_group, age,
    sbp,
    tcol_mmol,
    diab, 
    smoke)




#####Analysis: Globorisk Risk Prediction#####

#Apply the Function
globorisk_dummy <- globorisk_dummy %>%
  mutate(
    risk_recalib = globorisk_recalib_fast(
      sex = sex,
      age = age,
      sbp = sbp,
      tcol_mmol = tcol_mmol,
      diab = diab,
      smoke = smoke))


#Apply the Function
globorisk_dummy_mg_dl <- globorisk_dummy_mg_dl %>%
  mutate(
    risk_recalib = globorisk_recalib_fast(
      sex = sex,
      age = age,
      sbp = sbp,
      tcol_mmol = tcol_mmol,
      diab = diab,
      smoke = smoke))

#####Analysis: Table Creation (mmol/dl)#####

globorisk_pred <- globorisk_dummy %>%
  mutate(
    age_group = factor(
      age_group,
      levels = c("70-74", "65-69", "60-64", "55-59", "50-54", "45-49", "40-44"),
      ordered = TRUE),
    sbp = factor(
      sbp,
      levels = c(180, 160, 140, 120),
      ordered = TRUE),
    tcol_mmol = factor(
      tcol_mmol,
      levels = c(3, 4, 5, 6, 7),
      ordered = TRUE))


#Hombres
tab1.1<-apply_risk_colors(make_ft(make_globorisk_matrix(globorisk_pred, sex_value = 0, diab_value = 0,smoke_value =  0),
                          titulo = "Hombres | No diabéticas | No fumadoras"), 
                  make_globorisk_matrix(globorisk_pred, sex_value = 0, diab_value = 0,smoke_value =  0))

tab1.2<-apply_risk_colors(make_ft(make_globorisk_matrix(globorisk_pred, sex_value = 0, diab_value = 0,smoke_value =  1),
                                titulo = "Hombres | No diabéticas | Fumadoras"), 
                        make_globorisk_matrix(globorisk_pred, sex_value = 0, diab_value = 0,smoke_value =  1))

tab1.3<-apply_risk_colors(make_ft(make_globorisk_matrix(globorisk_pred, sex_value = 0, diab_value = 1,smoke_value =  0),
                                titulo = "Hombres | Diabéticas | No fumadoras"), 
                        make_globorisk_matrix(globorisk_pred, sex_value = 0, diab_value = 1,smoke_value =  0))

tab1.4<-apply_risk_colors(make_ft(make_globorisk_matrix(globorisk_pred, sex_value = 0, diab_value = 1,smoke_value =  1),
                                titulo = "Hombres | Diabéticas | Fumadoras"), 
                        make_globorisk_matrix(globorisk_pred, sex_value = 0, diab_value = 1,smoke_value =  1))


#Mujeres
tab2.1<-apply_risk_colors(make_ft(make_globorisk_matrix(globorisk_pred, sex_value = 1, diab_value = 0,smoke_value =  0),
                                  titulo = "Mujeres | No diabéticas | No fumadoras"), 
                          make_globorisk_matrix(globorisk_pred, sex_value = 1, diab_value = 0,smoke_value =  0))

tab2.2<-apply_risk_colors(make_ft(make_globorisk_matrix(globorisk_pred, sex_value = 1, diab_value = 0,smoke_value =  1),
                                  titulo = "Mujeres | No diabéticas | Fumadoras"), 
                          make_globorisk_matrix(globorisk_pred, sex_value = 1, diab_value = 0,smoke_value =  1))

tab2.3<-apply_risk_colors(make_ft(make_globorisk_matrix(globorisk_pred, sex_value = 1, diab_value = 1,smoke_value =  0),
                                  titulo = "Mujeres | Diabéticas | No fumadoras"), 
                          make_globorisk_matrix(globorisk_pred, sex_value = 1, diab_value = 1,smoke_value =  0))

tab2.4<-apply_risk_colors(make_ft(make_globorisk_matrix(globorisk_pred, sex_value = 1, diab_value = 1,smoke_value =  1),
                                  titulo = "Mujeres | Diabéticas | Fumadoras"), 
                          make_globorisk_matrix(globorisk_pred, sex_value = 1, diab_value = 1,smoke_value =  1))


tab1.1%>%flextable::save_as_docx(path="Tab_1.1.docx")
tab1.2%>%flextable::save_as_docx(path="Tab_1.2.docx")
tab1.3%>%flextable::save_as_docx(path="Tab_1.3.docx")
tab1.4%>%flextable::save_as_docx(path="Tab_1.4.docx")

tab2.1%>%flextable::save_as_docx(path="Tab_2.1.docx")
tab2.2%>%flextable::save_as_docx(path="Tab_2.2.docx")
tab2.3%>%flextable::save_as_docx(path="Tab_2.3.docx")
tab2.4%>%flextable::save_as_docx(path="Tab_2.4.docx")



#####Analysis: Table Creation (mg/dl)#####

globorisk_pred_mg_dl <- globorisk_dummy_mg_dl %>%
  mutate(
    age_group = factor(
      age_group,
      levels = c("70-74", "65-69", "60-64", "55-59", "50-54", "45-49", "40-44"),
      ordered = TRUE
    ),
    sbp = factor(
      sbp,
      levels = c(180, 160, 140, 120),
      ordered = TRUE
    ),
    tcol_mmol = factor(
      tcol_mmol,
      levels = c(2.6, 3.9, 5.2, 6.5, 7.8),
      labels = c("3", "4", "5", "6", "7"),
      ordered = TRUE
    ))


#Hombres
tab1.1_v2<-apply_risk_colors(make_ft(make_globorisk_matrix(globorisk_pred_mg_dl, sex_value = 0, diab_value = 0,smoke_value =  0),
                                  titulo = "Hombres | No diabéticas | No fumadoras"), 
                          make_globorisk_matrix(globorisk_pred_mg_dl, sex_value = 0, diab_value = 0,smoke_value =  0))

tab1.2_v2<-apply_risk_colors(make_ft(make_globorisk_matrix(globorisk_pred_mg_dl, sex_value = 0, diab_value = 0,smoke_value =  1),
                                  titulo = "Hombres | No diabéticas | Fumadoras"), 
                          make_globorisk_matrix(globorisk_pred_mg_dl, sex_value = 0, diab_value = 0,smoke_value =  1))

tab1.3_v2<-apply_risk_colors(make_ft(make_globorisk_matrix(globorisk_pred_mg_dl, sex_value = 0, diab_value = 1,smoke_value =  0),
                                  titulo = "Hombres | Diabéticas | No fumadoras"), 
                          make_globorisk_matrix(globorisk_pred_mg_dl, sex_value = 0, diab_value = 1,smoke_value =  0))

tab1.4_v2<-apply_risk_colors(make_ft(make_globorisk_matrix(globorisk_pred_mg_dl, sex_value = 0, diab_value = 1,smoke_value =  1),
                                  titulo = "Hombres | Diabéticas | Fumadoras"), 
                          make_globorisk_matrix(globorisk_pred_mg_dl, sex_value = 0, diab_value = 1,smoke_value =  1))


#Mujeres
tab2.1_v2<-apply_risk_colors(make_ft(make_globorisk_matrix(globorisk_pred_mg_dl, sex_value = 1, diab_value = 0,smoke_value =  0),
                                  titulo = "Mujeres | No diabéticas | No fumadoras"), 
                          make_globorisk_matrix(globorisk_pred_mg_dl, sex_value = 1, diab_value = 0,smoke_value =  0))

tab2.2_v2<-apply_risk_colors(make_ft(make_globorisk_matrix(globorisk_pred_mg_dl, sex_value = 1, diab_value = 0,smoke_value =  1),
                                  titulo = "Mujeres | No diabéticas | Fumadoras"), 
                          make_globorisk_matrix(globorisk_pred_mg_dl, sex_value = 1, diab_value = 0,smoke_value =  1))

tab2.3_v2<-apply_risk_colors(make_ft(make_globorisk_matrix(globorisk_pred_mg_dl, sex_value = 1, diab_value = 1,smoke_value =  0),
                                  titulo = "Mujeres | Diabéticas | No fumadoras"), 
                          make_globorisk_matrix(globorisk_pred_mg_dl, sex_value = 1, diab_value = 1,smoke_value =  0))

tab2.4_v2<-apply_risk_colors(make_ft(make_globorisk_matrix(globorisk_pred_mg_dl, sex_value = 1, diab_value = 1,smoke_value =  1),
                                  titulo = "Mujeres | Diabéticas | Fumadoras"), 
                          make_globorisk_matrix(globorisk_pred_mg_dl, sex_value = 1, diab_value = 1,smoke_value =  1))


tab1.1_v2%>%flextable::save_as_docx(path="Tab_V2_1.1.docx")
tab1.2_v2%>%flextable::save_as_docx(path="Tab_V2_1.2.docx")
tab1.3_v2%>%flextable::save_as_docx(path="Tab_V2_1.3.docx")
tab1.4_v2%>%flextable::save_as_docx(path="Tab_V2_1.4.docx")

tab2.1_v2%>%flextable::save_as_docx(path="Tab_V2_2.1.docx")
tab2.2_v2%>%flextable::save_as_docx(path="Tab_V2_2.2.docx")
tab2.3_v2%>%flextable::save_as_docx(path="Tab_V2_2.3.docx")
tab2.4_v2%>%flextable::save_as_docx(path="Tab_V2_2.4.docx")
