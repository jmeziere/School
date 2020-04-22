hist(RentClose$Price,
     xlab = "Rent Price Per Semester",
     ylab = "Number of Apartments",
     main = "Apartments by Price",
     col = c("skyblue","skyblue","skyblue","skyblue","gray","gray","gray","gray"))


plot_ly(data = RentClose, x = ~Price, y = ~WalkingMinutes, text = ~paste(Apartment)) %>%
  layout(title = "BYU-I Approved Housing",
         xaxis = list(title = "Price"),
         yaxis = list (title = "Walking Minutes"))

datatable(Rent, extensions="Responsive", options=list(lengthMenu=c(3,5,10)))