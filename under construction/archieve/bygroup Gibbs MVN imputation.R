data = obsdf
group = list_get_obs$list_obsdf[[pm]]$class

obj_mice = mice(data, m = 5, method = "bygroup", group = group, 
                imputationFunction = "norm")

obj_mice = mice(data, m=5, blocks = list(block1 = c("Y1","Y2","Y3","Y4","Y5","Xinc1", "Xinc2")), 
                method = "bygroup",group = group, 
                data.init =list_get_comp$dfcom, maxit = 100)


obj_mice = mice(data, m=5, method = "bygroup",group = group,imputationFunction = "norm.boot", 
                data.init =list_get_comp$dfcom, maxit = 20)
