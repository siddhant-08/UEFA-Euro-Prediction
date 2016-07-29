library(data.table)
library(dplyr)

df<- fread('matches.csv', header = T)
df<- as.data.frame(df)
# View(df)
teams<- c('Albania','Austria','Belgium','Croatia','Czech Republic','England',
          'France','Germany','Hungary','Iceland','Italy','North. Ireland',
          'Poland','Portugal','Ireland','Romania','Russia','Slovakia','Spain',
          'Sweden','Switzerland','Turkey','Ukraine','Wales')

# all matches between relevant teams
# 771
euro_df<- filter(df, home_team %in% teams, away_team %in% teams)

euro_df$score<-''

do_score<- function(home,away){
  return(paste0(home,'-',away))
}

euro_df$score<- mapply(do_score, euro_df$home_goals, euro_df$away_goals)

euro_df$score<- as.factor(euro_df$score)
#
# for checking distribution. the original model slightly inflated likelihood 
# for 0-0, 0-1 and 1-0 scores.
barplot(table(euro_df$score))
#,'QE-12','QE-13','QE-14'

valid_matches<- c('EC-12','EC-13','EC-14','QE-12','QE-13','QE-14','WC-20')
euro_df<- filter(euro_df, type %in% valid_matches)

euro_df$home_team<- as.factor(euro_df$home_team)
euro_df$away_team<- as.factor(euro_df$away_team)


euro_df$months<-0

find_months<- function(year,month,ref_year,ref_month){
  return( 12*(year - ref_year) + (month - ref_month))
}

euro_df$months<-mapply(find_months,euro_df$Year, euro_df$Month,2002,1)

euro_df$months<- sapply(euro_df$months, function(x,max_month){ x<- max_month - x + 1}, max_month=euro_df$months[nrow(euro_df)] )


#valid_matches<- c('EC-12','EC-13','EC-14','QE-12','QE-13','QE-14')
#prof_euro_df<- filter(euro_df, type %in% valid_matches)
#euro_df<- prof_euro_df

likelihood <- function(y1, y2, lambda, mu, months,phi=0,rho=0){
  #rho=0, independence
  #y1: home goals
  #y2: away goals=
  #exp(-phi*(months) )
  sum( exp(-phi*months) * (log(dpois(y1, lambda) ) + log(dpois(y2, mu) ) ) ) 
}

create_base <- function(euro_df){
  
  hm <- model.matrix(~ home_team - 1, data=euro_df,contrasts.arg=list(home_team='contr.treatment'))
  am <- model.matrix(~ away_team -1, data=euro_df)
  
  team.names <- unique(c(levels(euro_df$home_team), levels(euro_df$away_team)))
  
  return(list(
    homeTeamDM=hm,
    awayTeamDM=am,
    homeGoals=euro_df$home_goals,
    awayGoals=euro_df$away_goals,
    teams=team.names
  )) 
}

DCoptimFn <- function(params, DCm, months,phi.p=0){
  
  home.p <- params[1]
  
  nteams <- length(DCm$teams)
  attack.p <- matrix(params[2 : ( nteams + 1)], ncol=1)
  defence.p <- matrix(params[( nteams + 2) : length(params) ], ncol=1)

  lambda <- exp(DCm$homeTeamDM %*% attack.p + DCm$awayTeamDM %*% defence.p + home.p)
  mu <- exp(DCm$awayTeamDM %*% attack.p + DCm$homeTeamDM %*% defence.p)
  
  return(
    likelihood(DCm$homeGoals, DCm$awayGoals, lambda, mu, months, phi.p, rho.p) * -1
  )
}

attack_constraint <- function(params, DCm, ...){
  nteams <- length(DCm$teams)
  attack.p <- matrix(params[2 : ( nteams+1 )], ncol=1)
  return((sum(attack.p) / nteams) - 1)
}

dcm<- create_base(euro_df)
#initial parameter estimates
attack_params <- rep(.5, times=nlevels(as.factor(euro_df$home_team)))
defence_params <- rep(-0.25, times=nlevels(as.factor(euro_df$home_team)))
home_param <- 0.2
phi <-0.03
init_params <- c(home_param,attack_params, defence_params)
names(init_params) <- c('HOME', paste('Attack', dcm$teams, sep='.'), paste('Defence', dcm$teams, sep='.'))


library(alabama)
res <- auglag(par=par.inits, 
              fn=DCoptimFn,
              heq=attack_constraint,
              DCm=dcm, 
              months= euro_df$months,
              phi.p=phi)

match<- fread('predict.csv',header = T)
fifa<- fread('fifa.csv',header=T)
match<- as.data.frame(match)

matches<- match
# adding columns to the table
matches$lambda<-0
matches$mu<-0


calc_lambda<- function(home_team,away_team) {
  home_data<-fifa[fifa$Team==home_team]
  away_data<-fifa[fifa$Team==away_team]
  
  lambda<- 0.7*(res$par[paste0('Attack.',home_team)] + res$par[paste0('Defence.',away_team)]) + 
    0.2*((home_data$rating - away_data$rating)/717) + 
    0.1*((0.6*home_data$GF + 0.4*away_data$GA)/25.4)
  # 1384 - 667
  # 33*0.6 + 14*0.4
  
  mu<- 0.7*(res$par[paste0('Attack.',away_team)] + res$par[paste0('Defence.',home_team)]) + 
    0.2*((away_data$rating - home_data$rating)/717) + 
    0.1*((0.6*away_data$GF + 0.4*home_data$GA)/25.4)
  
  if(home_team =='France' || home_team=='Belgium'){
    lambda<- lambda + res$par['HOME']
  }
  
  if(away_team =='France' || away_team=='Belgium'){    # to account for Belgium's stellar rise in recent years
    mu<- mu + res$par['HOME']
  }
  
  return (list(exp(lambda), exp(mu)))
}

la_mu<-mapply(calc_lambda, matches$home_team,matches$away_team)

matches$lambda<- la_mu[1,]
matches$mu <-la_mu[2,]


matches$home_score<-0.0
matches$away_score<-0.0

find_score<- function(lambda, mu,group_stage){
  maxgoal <- 8 # will be useful later
  probability_matrix <- dpois(0:maxgoal, lambda) %*% t(dpois(0:maxgoal, mu))
  
  # Predicting the most probable non-draw score for knockouts.
  # Not enough data to predict penalty scores
  
  if(group_stage==0){
    diag(probability_matrix)<-0
  }
  
  score<-as.data.frame( which( probability_matrix==max(probability_matrix), arr.ind = T ))
  return (score)
}

find_odds<- function(lambda, mu){
  maxgoal <- 8 
  probability_matrix <- dpois(0:maxgoal, lambda) %*% t(dpois(0:maxgoal, mu))
  
  draw<-0
  win<-0
  loss<-0
  win <- sum(probability_matrix[lower.tri(probability_matrix)])
  draw <- sum(diag(probability_matrix))
  loss <- sum(probability_matrix[upper.tri(probability_matrix)])
  return (c(win,draw,loss))
}


score<- mapply(find_score, matches$lambda, matches$mu,matches$group_stage)
matches$home_score<-as.numeric(score[1,])
matches$away_score<- as.numeric(score[2,])
matches$home_score<- sapply(matches$home_score,FUN = function(x){x= x-1 })
matches$away_score<- sapply(matches$away_score,FUN = function(x){x= x-1 })

odds<- mapply(find_odds,matches$lambda, matches$mu)
matches$win<- odds[1,]
matches$draw<-odds[2,]
matches$loss<- odds[3,]

View(matches)
