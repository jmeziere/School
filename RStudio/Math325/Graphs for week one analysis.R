plot_ly(data = RentClose, x = ~Price, y = ~WalkingMinutes, text = ~paste(Apartment)) %>%
  layout(title = "BYU-I Approved Housing",
         xaxis = list(title = "Price"),
         yaxis = list (title = "Walking Minutes"))
datatable(Rent, extensions="Responsive", options=list(lengthMenu=c(3,5,10)))
