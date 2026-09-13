# Record the R environment

sink(file.path("environment", "sessionInfo.txt"))
print(sessionInfo())
sink()
