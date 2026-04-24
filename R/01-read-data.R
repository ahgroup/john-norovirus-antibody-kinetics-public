###
# 01-read-data.R
###

func_process_ni_data <- function(data_path) {
  
  data <- readRDS(data_path)
  
  return(data)
}

func_process_nv_data <- function(data_path) {
  
  data <- readRDS(data_path)
  
  return(data)
}

# supplement ====
# change lower censoring bound
func_alt_process_ni_data <- function(data){
  
  # increase lower bound
  dat_increase <- data %>%
    mutate(
      y = case_when(
        cens == 2 ~ 10e-3, # 10^(-2)
        TRUE ~ y
      ),
      logy = log(y)
    )
  
  # decrease lower bound
  dat_decrease <- data %>%
    mutate(
      y = case_when(
        cens == 2 ~ 10e-15, # 10^(-14)
        TRUE ~ y
      ),
      logy = log(y)
    )
  
  list(
    dat_ni_increase = dat_increase,
    dat_ni_decrease = dat_decrease
  )
  
}

func_alt_process_nv_data <- function(data){
  
  dat_increase <- data %>%
    mutate(
      y = case_when(
        cens == 2 ~ 10e-3,
        TRUE ~ y
      ),
      logy = log(y)
    )
  
  # decrease lower bound
  dat_decrease <- data %>%
    mutate(
      y = case_when(
        cens == 2 ~ 10e-15,
        TRUE ~ y
      ),
      logy = log(y)
    )
  
  list(
    dat_nv_increase = dat_increase,
    dat_nv_decrease = dat_decrease
  )
  
}

# END OF SCRIPT ====