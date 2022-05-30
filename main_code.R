library(sfcr)
library(tidyverse)


### Calibragem
calib_sm_inv <- sfcr_set(
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
  S ~ (Y*(1+g))/(1+iota-beta),
  INV ~ beta*S,
  Fb ~ i*((L+INV-D)/(1+g)),
  #Firmas
  W ~ (1 - pi) * Y,
  I ~ h*Y,
  delta ~ 1 - (K-I)*(1+g) / K,
  Ft ~ S - (1-pi)*Y + (INV - INV/(1+g)) - i*L/(1+g) - i*INV/(1+g),
  Fu ~ sf*Ft,
  Fd ~ Ft-Fu,
  pi ~ mu/(1+mu),
  G ~ sigma * S*((1+g)^-1),
  Tg ~ tau * (W + Fd +i*((B+D)/(1+g))),
  C ~ S - G - I,
  alpha2 ~ (C - alpha1*(1- tau)*W) * (1+g)/Vh,
  r_n ~ Ft / (K/(1+g)),
  r_g ~ S/(K/(1+g)), 
)

sm_inv_vals <- sfcr_set(
  K ~ 100,
  tau ~ 0.37,
  sigma ~ 0.34,
  alpha1 ~ 0.8,
  u ~ 0.8,
  a ~ 0.1,
  x ~ 0.001,
  lambda0 ~ 0.08,
  v ~ 2.5,
  i ~ 0.02,
  mu ~ 0.7,
  sf ~ 0.4,
  h ~ 0.2,
  g ~ 0.02,
  l ~ 0.80,
  b ~ 0.75,
  savf_ratio ~ -0.02,
  beta ~ 0.2,
  iota ~ 0.2,
)

sm_inv_ini_vals <- sfcr_baseline(
  equations = calib_sm_inv,
  external = sm_inv_vals,
  periods = 2) %>%
  filter(row_number() == 2)

#############

f_d_h <- function(u, un, x, h, gh,gamma){
  if (abs(u-un) > x){
    (h/(1+gh))*gamma*(u-un)
  }else {
    0
  }
}

model_sm_inv <-sfcr_set(
  # Total da Economia
  S ~ C + I + G,
  
  #Bancos
  Fb ~ i*L[-1] + i*INV[-1] - i*D[-1],
  
  # Governo
  B ~ B[-1] + G - Tg + i*B[-1],
  G ~ sigma * Y[-1],
  Tg ~ tau * Y_h,
  # Famílias
  Y_h ~ W + Fd + i*(B[-1] + D[-1]) + Fb,
  W ~ (1 - pi) * Y,
  Y_d ~ (1-tau) * Y_h,
  C ~ alpha1*(1- tau)*W + alpha2*Vh[-1],
  S_h ~ Y_d - C,
  lambda ~ lambda0 - i, 
  pe ~ lambda * Vh / E,
  Vh ~ Vh[-1] + S_h + (pe - pe[-1]) * E[-1],
  D ~ D[-1] + S_h - pe * (E - E[-1]) - (B - B[-1]),
  #Firmas
  #Equações dos Estoques
  Y ~ (1+iota)*S_e - INV[-1],
  S_e ~ S[-1],
  INV ~ Y - S + INV[-1],
  ######
  pi ~ mu / (1+mu),
  gh ~ (h - h[-1])/h[-1],
  d_h ~ f_d_h(u, un, x, h, gh,gamma),
  h ~ h[-1] + d_h,
  I ~ h*Y,
  K ~ K[-1] - delta*K[-1] + I,
  Y_fc ~ K[-1]/v,
  u ~ Y / Y_fc,
  g_k ~ (h*u /v) - delta,
  E ~ a*K[-1],
  L ~ L[-1] + I - Fu - pe*(E - E[-1]) + (INV - INV[-1]),
  Ft ~ S - (1-pi)*Y - i*L[-1] + (INV - INV[-1]) - i*INV,
  Fu ~ sf*Ft,
  Fd ~ (1-sf)*Ft,
  F_g ~ pi*Y,
  r_n ~ Ft/K[-1],
  r_g ~ S/K[-1], 
######################  
g_y ~ -1 + Y/Y[-1],
g_inv ~ -1 + INV/INV[-1],
g_vh ~ -1 + Vh/Vh[-1],
g_sale ~ -1 + S/S[-1],
l ~ L/K[-1],
b ~ B/K[-1],
vh ~ Vh/K[-1],
eq_value ~ pe*E/K[-1],
z ~ alpha2*Vh[-1]/Y
)

sm_inv_pars <- sfcr_set(
  tau ~ sm_inv_ini_vals$tau,
  sigma ~ sm_inv_ini_vals$sigma,
  alpha1 ~ sm_inv_ini_vals$alpha1,
  alpha2 ~ sm_inv_ini_vals$alpha2,
  a ~ sm_inv_ini_vals$a,
  iota ~ sm_inv_ini_vals$iota,
  x ~ 0.001,
  lambda0 ~ sm_inv_ini_vals$lambda0,
  gamma ~ 0.014,
  delta ~ sm_inv_ini_vals$delta,
  i ~ sm_inv_ini_vals$i,
  mu ~ sm_inv_ini_vals$mu,
  sf ~ 0.4,
  un ~ sm_inv_ini_vals$u,
  v ~ 2.5,
)

sm_inv_ini <- sfcr_set(
  Y ~ sm_inv_ini_vals$Y,
  B ~ sm_inv_ini_vals$B,
  D ~ sm_inv_ini_vals$D,
  Vh ~ sm_inv_ini_vals$Vh,
  h ~ sm_inv_ini_vals$h,
  K ~ sm_inv_ini_vals$K,
  L ~ sm_inv_ini_vals$L,
  E ~ sm_inv_ini_vals$E,
  pe ~ sm_inv_ini_vals$pe,
  u ~ sm_inv_ini_vals$u,
  S ~ sm_inv_ini_vals$S,
  INV ~ sm_inv_ini_vals$INV,
)

sm_inv_baseline <- sfcr_baseline( 
  equations = model_sm_inv,
  external = sm_inv_pars,
  initial = sm_inv_ini,
  periods = 500)

####
a <- sm_inv_baseline[nrow(sm_inv_baseline),] %>% select(-period)
a <- a %>% slice(rep(1:n(), each = as.numeric(nrow(sm_inv_baseline))))
a$period <- 1:as.numeric(nrow(sm_inv_baseline))
sm_inv_base_comp <- a
sm_inv_base_comp$Scenario <- "Baseline"
####


#########################
### Modificação em iota #
#########################

iota_schock_sm_inv_baseline <- sfcr_shock(
  variables = sfcr_set(
    iota ~ sm_inv_baseline$iota[2] + 0.03),
  start = 1, 
  end = 500)

sm_inv_ext_1 <- sfcr_scenario(
  baseline = sm_inv_baseline,
  scenario = iota_schock_sm_inv_baseline,
  periods = 500)

reshaped_sm_inv_ext_1 <- sm_inv_ext_1 %>%
  pivot_longer(cols = -period )

sm_inv_ext_1$Scenario <- "Choque em iota"
sm_inv_scenario_1 <- rbind(sm_inv_base_comp,
                           sm_inv_ext_1)

###∆ vs baseline

a <-sm_inv_scenario_1 %>% 
  ggplot( aes( x = period, y = u , group = Scenario, 
               color = Scenario, linetype = Scenario)) +geom_line( ) +
  theme(legend.position = "bottom", legend.title = element_blank()) + xlab( "Tempo" )

b <- sm_inv_scenario_1 %>% 
  ggplot( aes( x = period, y = h, group = Scenario, 
               color = Scenario, linetype = Scenario)) +geom_line( ) +
  theme(legend.position = "bottom", legend.title = element_blank()) + xlab( "Tempo" )

cowplot::plot_grid(a,b, labels = c("A","B"))

c <- sm_inv_scenario_1 %>% 
  ggplot( aes( x = period, y = vh, group = Scenario, 
               color = Scenario, linetype = Scenario)) +geom_line( ) +
  theme(legend.position = "bottom", legend.title = element_blank()) + xlab( "Tempo" )

d <- sm_inv_scenario_1 %>% 
  ggplot( aes( x = period, y = l, group = Scenario, 
               color = Scenario, linetype = Scenario)) +geom_line( )+
  theme(legend.position = "bottom", legend.title = element_blank()) + xlab( "Tempo" )


e <- sm_inv_scenario_1 %>% 
  ggplot( aes( x = period, y = b, group = Scenario, 
               color = Scenario, linetype = Scenario)) +geom_line( )+
  theme(legend.position = "bottom", legend.title = element_blank()) + xlab( "Tempo" )

f1 <-sm_inv_ext_1 %>% select("eq_value","period")
f_eq <- as.data.frame(f1[,1]/as.numeric(f1[1,1]))
f_eq$period <- 1:nrow(sm_inv_ext_1)

f <- f_eq %>%
  ggplot( aes( x = period, y = eq_value)) +
  geom_line( ) +
  xlab( "Tempo" )

cowplot::plot_grid(c,d,e,f, labels = c("A","B","C","D"))

###∆
a <- reshaped_sm_inv_ext_1 %>% 
  filter( name %in% c( "g_sale", "g_y") ) %>%
  ggplot( aes( x = period, y = value, group = name, linetype = name, color = name ) ) +
  geom_line() +
  labs(y=" ", x= "Tempo") +
  theme(legend.position = "bottom", legend.title = element_blank()) + xlab( "Tempo" )

b <- reshaped_sm_inv_ext_1 %>% 
  filter( name %in% c( "g_k", "g_vh") ) %>%
  ggplot( aes( x = period, y = value, group = name, linetype = name, color = name ) ) +
  geom_line() +
  labs(y=" ", x= "Tempo") +
  theme(legend.position = "bottom", legend.title = element_blank()) + xlab( "Tempo" )

cowplot::plot_grid(a,b, labels = c("A","B"))

### variáveis normalizadas

f1 <-sm_inv_ext_1 %>% select("r_g","r_n","period")
r_g <- f1[,1]/with(sm_inv_baseline[nrow(sm_inv_baseline),],r_g)
r_n <- f1[,2]/with(sm_inv_baseline[nrow(sm_inv_baseline),],r_n)
r_norm <- cbind(r_g,r_n)
r_norm$period <- f1$period

b <-r_norm %>% pivot_longer(cols = -period ) %>%
  filter( name %in% c( "r_n", "r_g") ) %>%
  ggplot( aes( x = period, y = value, group = name, linetype = name, color = name ) ) +
  geom_line() +
  labs(y=" ", x= "Tempo") +
  theme(legend.position = "bottom", legend.title = element_blank()) + xlab( "Tempo" )

plot(b)

##########
#Gráficos do Baseline
##########

reshaped_sm_inv <- sm_inv_baseline %>%
  pivot_longer(cols = -period )

a <-reshaped_sm_inv %>% 
  filter( name %in% c( "Y", "INV") ) %>% 
  ggplot( aes( x = period, y = value, group = name, linetype = name, color = name ) ) +
  geom_line() +
  labs(y=" ", x= "Tempo") +
  theme(legend.position = "bottom", legend.title = element_blank()) 

b <-reshaped_sm_inv %>% 
  filter( name %in% c( "g_y", "g_inv") ) %>% 
  ggplot( aes( x = period, y = value, group = name, linetype = name, color = name ) ) +
  geom_line() +
  labs(y=" ", x= "Tempo") +
  theme(legend.position = "bottom", legend.title = element_blank()) 

cowplot::plot_grid(a,b, labels = c("A","B"))

a<-sm_inv_baseline %>%
  ggplot( aes( x = period, y = INV/S) ) +
  geom_line( ) +
  xlab( "Tempo" )

b<-sm_inv_baseline %>%
  ggplot( aes( x = period, y = (INV-lag(INV))/S) ) +
  geom_line( ) +
  xlab( "Tempo" ) +
  ylab("d_INV/Y")

cowplot::plot_grid(a,b, labels = c("A","B"))


####################################################
### Modificação na Propensão a Consumir da Riqueza #
####################################################

alpha2_mod_sm_inv_baseline <- sfcr_shock(
  variables = sfcr_set(
    alpha2 ~ sm_inv_baseline$alpha2[2] + 0.002),
  start = 1, 
  end = 500)

sm_inv_ext_2 <- sfcr_scenario(
  baseline = sm_inv_baseline,
  scenario = alpha2_mod_sm_inv_baseline,
  periods = 500)

reshaped_sm_inv_ext_2 <- sm_inv_ext_2 %>%
  pivot_longer(cols = -period )

sm_inv_ext_2$Scenario <- "Choque em alpha_2"
sm_inv_scenario_2 <- rbind(sm_inv_base_comp,
                           sm_inv_ext_2)

###∆ vs baseline

a <-sm_inv_scenario_2 %>% 
  ggplot( aes( x = period, y = u , group = Scenario, 
               color = Scenario, linetype = Scenario)) +geom_line( ) +
  theme(legend.position = "bottom", legend.title = element_blank()) + xlab( "Tempo" )

b <- sm_inv_scenario_2 %>% 
  ggplot( aes( x = period, y = h, group = Scenario, 
               color = Scenario, linetype = Scenario)) +geom_line( ) +
  theme(legend.position = "bottom", legend.title = element_blank()) + xlab( "Tempo" )

c <-sm_inv_scenario_2 %>% 
  ggplot( aes( x = period, y = z, group = Scenario, 
               color = Scenario, linetype = Scenario)) +geom_line( ) +
  theme(legend.position = "bottom", legend.title = element_blank()) + xlab( "Tempo" )

d <- reshaped_sm_inv_ext_2 %>% 
  filter( name %in% c( "g_k", "g_vh") ) %>%
  ggplot( aes( x = period, y = value, group = name, linetype = name, color = name ) ) +
  geom_line() +
  labs(y=" ", x= "Tempo") +
  theme(legend.position = "bottom", legend.title = element_blank()) + xlab( "Tempo" ) 

cowplot::plot_grid(a,b,c,d, labels = c("A","B","C","D"))



a <-sm_inv_scenario_2 %>% 
  ggplot( aes( x = period, y = l , group = Scenario, 
               color = Scenario, linetype = Scenario)) +geom_line( ) +
  theme(legend.position = "bottom", legend.title = element_blank()) + xlab( "Tempo" )

b <- sm_inv_scenario_2 %>% 
  ggplot( aes( x = period, y = b, group = Scenario, 
               color = Scenario, linetype = Scenario)) +geom_line( ) +
  theme(legend.position = "bottom", legend.title = element_blank()) + xlab( "Tempo" )

c <-sm_inv_scenario_2 %>% 
  ggplot( aes( x = period, y = vh, group = Scenario, 
               color = Scenario, linetype = Scenario)) +geom_line( ) +
  theme(legend.position = "bottom", legend.title = element_blank()) + xlab( "Tempo" )

f2 <-sm_inv_ext_2 %>% select("eq_value","period")
f_eq <- as.data.frame(f2[,1]/as.numeric(f2[1,1]))
f_eq$period <- 1:nrow(sm_inv_ext_1)

d  <- f_eq %>%
  filter(period <= 50) %>%
  ggplot( aes( x = period, y = eq_value)) +
  geom_line( ) +
  xlab( "Tempo" )

cowplot::plot_grid(a,b,c,d, labels = c("A","B","C","D"))

f2 <-sm_inv_ext_2 %>% select("r_g","r_n","period")
r_g <- f2[,1]/with(sm_inv_baseline[nrow(sm_inv_baseline),],r_g)
r_n <- f2[,2]/with(sm_inv_baseline[nrow(sm_inv_baseline),],r_n)
r_norm_2 <- cbind(r_g,r_n)
r_norm_2$period <- 1:nrow(sm_inv_ext_2)


a<-r_norm_2 %>% pivot_longer(cols = -period ) %>%
  filter( name %in% c( "r_n", "r_g") ) %>% 
  ggplot( aes( x = period, y = value, group = name, linetype = name, color = name ) ) +
  geom_line() +
  labs(y=" ", x= "Tempo") +
  theme(legend.position = "bottom", legend.title = element_blank()) + xlab( "Tempo" ) 



b<-reshaped_sm_inv_ext_2 %>% 
  filter( name %in% c( "g_sale", "g_y") ) %>%
  ggplot( aes( x = period, y = value, group = name, linetype = name, color = name ) ) +
  geom_line() +
  labs(y=" ", x= "Tempo") +
  theme(legend.position = "bottom", legend.title = element_blank()) + xlab( "Tempo" )

cowplot::plot_grid(a,b, labels = c("A","B"))