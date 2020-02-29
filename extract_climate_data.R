# Extract climate data from NOAA product

library("tabulizer")

Climdata <- extract_tables("Data/Climate_data_West_Thompson_Lake.pdf", pages = 1, 
                           area = list(c(10,149,100,200))); View(Climdata[[1]])

Climdata <- extract_areas("Data/Climate_data_West_Thompson_Lake.pdf", pages = 1)


get_page_dims("Data/Climate_data_West_Thompson_Lake.pdf", pages = 1)
