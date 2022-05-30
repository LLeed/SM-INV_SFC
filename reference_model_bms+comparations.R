library(sfcr)
library(tidyverse)


### Calibragem 
calib_lmc <- sfcr_set(
  L ~ l * K,
  B ~ b * K,
  D ~ L,
  SAVg ~ -g * B / ( 1 + g ),
  SAVf ~ savf_ratio * K,
  SAVh ~ -SAVg - SAVf,
  Vh ~ SAVh * ( 1 + g ) / g,
  lambda ~ lambda0 - i,
  E ~ a * (K*((1+g)^-1)),
  pe ~ lambda * Vh / E,
  Y ~ (u*K)/(v*(1+g)),
  W ~ (1 - pi) * Y,
  I ~ h*Y, 
  delta ~ 1 - (K-I)*(1+g) / K,
  Ft ~ pi*Y - i*L/(1+g),
  Fu ~ sf*Ft,
  Fd ~ Ft-Fu,pi ~ mu/(1+mu),
  G ~ sigma * Y*((1+g)^-1),
  Tg ~ tau * (W + Fd +i*((B+D)/(1+g))),
  C ~ Y - G - I,
  alpha2 ~ (C - alpha1*(1- tau)*W) * (1+g)/Vh,
)

bms_vals <- sfcr_set(
  K ~ 100,
  tau ~ 0.37,
  sigma ~ 0.34,
  alpha1 ~ 0.8,
  u ~ 0.8,
  a ~ 0.1,
  x ~ 0.001,
  lambda0 ~0.08,
  v ~ 2.5,
  i ~ 0.02,
  mu ~ 0.7,
  sf ~ 0.4,
  h ~ 0.2,
  g ~ 0.02,
  l ~ 0.80,
  b ~ 0.75,
  savf_ratio ~ -0.02,
)


bms_ini_vals <- sfcr_baseline(
  equations = calib_lmc,
  external = bms_vals,
  periods = 2, method = "Gauss") %>%
  filter(row_number() == 2)

f_d_h <- function(u, un, x, h, gh,gamma){
  if (abs(u-un) > x){ #x é o valor do corredor
    (h/(1+gh))*gamma*(u-un)
  }else {
    0
  }
}

model_lmc <-sfcr_set(
  Y ~ C + I + G,
  B ~ B[-1] + G - Tg + i*B[-1],
  G ~ sigma * Y[-1],
  Tg ~ tau * Y_h,
  Y_h ~ W + Fd + i*(B[-1] + D[-1]),
  W ~ (1 - pi) * Y,
  Y_d ~ (1-tau) * Y_h,
  C ~ alpha1*(1- tau)*W + alpha2*Vh[-1],
  S_h ~ Y_d - C,
  lambda ~ lambda0 - i, 
  pe ~ lambda * Vh / E,
  Vh ~ Vh[-1] + S_h + (pe - pe[-1]) * E[-1],
  D ~ D[-1] + S_h - pe * (E - E[-1]) - (B - B[-1]),
  #Firmas
  pi ~ mu / (1+mu),
  gh ~ (h - h[-1])/h[-1],
  d_h ~ f_d_h(u, un, x, h, gh,gamma),
  h ~ h[-1] + d_h,
  I ~ h*Y,
  K ~ K[-1] - delta*K[-1] + I,
  Y_fc ~ K[-1]/v,
  u ~ Y / Y_fc,
  g_k ~ (h*u /v) - delta,
  L ~ L[-1] + I - Fu - pe*(E - E[-1]),
  E ~ a*K[-1],
  Ft ~ pi*Y - i*L[-1],
  Fu ~ sf*Ft,
  Fd ~ (1-sf)*Ft,
  F_g ~ pi*Y,
  r_n ~ (pi* u / v) - i*L[-1]/K[-1],
  r_g ~ pi*u/v
)

bms_pars <- sfcr_set(
  tau ~ bms_ini_vals$tau,
  sigma ~ bms_ini_vals$sigma,
  alpha1 ~ bms_ini_vals$alpha1,
  alpha2 ~ bms_ini_vals$alpha2,
  a ~ bms_ini_vals$a,
  x ~ 0.001,
  lambda0 ~ bms_ini_vals$lambda0,
  gamma ~ 0.014,
  delta ~ bms_ini_vals$delta,
  i ~ bms_ini_vals$i,
  v ~ 2.5,
  mu ~ bms_ini_vals$mu,
  sf ~ 0.4,
  un ~ bms_ini_vals$u,
)

bms_ini <- sfcr_set(
  Y ~ bms_ini_vals$Y,
  B ~ bms_ini_vals$B,
  D ~ bms_ini_vals$D,
  Vh ~ bms_ini_vals$Vh,
  h ~ bms_ini_vals$h,
  K ~ bms_ini_vals$K,
  L ~ bms_ini_vals$L,
  E ~ bms_ini_vals$E,
  pe ~ bms_ini_vals$pe,
  u ~ bms_ini_vals$u,
)


bms_baseline <- sfcr_baseline( 
  equations = model_lmc,
  external = bms_pars,
  initial = bms_ini,
  periods = 500)

########
#Gráficos para comparação - OBS: O baseline do modelo principal deve ser
#simulado antes
########

a_1 <- bms_baseline %>%
  select("period","Y", "u", "g_k")

a_2 <- sm_inv_baseline %>%
  select("period","Y","u","g_k")

a_1$Scenario <- "Brochier e Silva (2019)"
a_2$Scenario <- "Baseline"

sim_scenario <- rbind(a_1, a_2)

sim_scenario %>%
  ggplot( aes( x = period, y = u/0.8  , group = Scenario, 
               color = Scenario, linetype = Scenario))+
  geom_line( ) + labs(y="u/un", x= "Tempo") + xlab("Tempo") +
  theme(legend.position = "bottom", legend.title = element_blank())
