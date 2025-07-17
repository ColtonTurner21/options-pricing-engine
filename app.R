# Black-Scholes Formula for European options
black_scholes <- function(S, K, r, sigma, T, type = c("call", "put")) {
  type <- match.arg(type)
  d_1 <- (log(S/K) + (r + 0.5 * sigma^2) * T) / (sigma * sqrt(T))
  d_2 <- d_1 - sigma * sqrt(T)
  if (type == "call") {
    price <- S * pnorm(d_1) - K * exp(-r * T) * pnorm(d_2)
  } else { 
    price <- -K * exp(-r * T) * pnorm(-d_2) - S * pnorm(-d_1) 
  }
  return(price)
}
black_scholes(100, 100, 0.05, 0.2, 1, "call")

# Binomial Tree Method for American options
binomial_american <- function(S, K, r, sigma, T, N = 100, type = c("call", "put")) {
  type <- match.arg(type)
  dt <- T/N
  u <- exp(sigma * sqrt(dt))
  d <- 1/u
  p <- (exp(r * dt) - d) / (u-d)
  discount <- exp(-r * dt)
  
  # Option values at maturity
  asset_prices <- S * u^(0:N) * d^((N):0)
  if (type == "call")  {
    values <- pmax(asset_prices - K, 0)
  } else {
    values <- pmax(K - asset_prices, 0)
  }
  
  # Backwards induction
  for (i in (N-1) : 0) {
    asset_prices <- S * u^(0:i) * d^((i):0)
    values <- discount * (p * values[2:(i+2)] + (1-p) * values[1:(i+1)])
    if (type == "call")  {
      values <- pmax(values, asset_prices - K)
    } else {
      values <- pmax(values, K - asset_prices)
    }
  }
  return(values[1])
}

# Monte Carlo Simulation for Asian options
asian_option_monte_carlo <- function(S, K, r, sigma, T, n_sim = 10000, 
                                     n_steps = 100, type = c("call", "put"))  {
  type <- match.arg(type)
  dt <- T/n_steps
  discount <- exp(-r * T)
  payoffs <- numeric(n_sim)
  for (i in 1:n_sim) {
    prices <- numeric(n_steps + 1)
    prices[1] <- S
    for (j in 2:(n_steps + 1)) {
      z <- rnorm(1)
      prices[j] <- prices[j-1] * exp((r - 0.5 * sigma^2) * dt + sigma * sqrt(dt)
                                     * z) }
    average_price <- mean(prices)
    if (type == "call") {
      payoffs[i] <- max(average_price - K, 0)
    } else {
      payoffs[i] <- max(K - average_price, 0)
    }
  }
  price <-discount * mean(payoffs)
  return(price)
}
# Shiny Dashboard 

library(shiny)
ui <- fluidPage(
  titlePanel("Options Pricing Dashboard"),
  sidebarLayout(
    sidebarPanel(
      selectInput("optionType", "Option Type", 
                  choices = c("European Call", "European Put", 
                              "American Call", "American Put",
                              "Asian Call", "Asian Put")),
      numericInput("S", "Stock Price (S)", value = 100, min = 1),
      numericInput("K", "Strike Price (K)", value = 100, min = 1),
      numericInput("r", "Risk-Free Rate (r)", value = 0.05, min = 0, step = 0.01),
      numericInput("sigma", "Volatility (sigma)", value = 0.2, min = 0.01, step = 0.01),
      numericInput("T", "Time to Expiration (Years)", value = 1, min = 0.01, step = 0.01),
      actionButton("goButton", "Price Option")
    ),
    mainPanel(
      h3("Option Price:"),
      verbatimTextOutput("optionPrice")
    )
  )
)

server <- function(input, output) {
  observeEvent(input$goButton, {
    output$optionPrice <- renderPrint({
      S <- input$S; K <- input$K; r <- input$r; sigma <- input$sigma; T <- input$T
      type <- input$optionType
      if (type == "European Call") {
        price <- black_scholes(S, K, r, sigma, T, type = "call")
      } else if (type == "European Put") {
        price <- black_scholes(S, K, r, sigma, T, type = "put")
      } else if (type == "American Call") {
        price <- binomial_american(S, K, r, sigma, T, type = "call")
      } else if (type == "American Put") {
        price <- binomial_american(S, K, r, sigma, T, type = "put")
      } else if (type == "Asian Call") {
        price <- asian_option_monte_carlo(S, K, r, sigma, T, type = "call")
      } else if (type == "Asian Put") {
        price <- asian_option_monte_carlo(S, K, r, sigma, T, type = "put")
      }
      price
    })
  })
}

shinyApp(ui = ui, server = server)