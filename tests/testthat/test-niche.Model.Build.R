# Test Start ###
test_that("format testing on the built model1",{
  data(en.vir)
  data(bak.vir)
  data<-data.frame(species=rep("Acosmeryx anceus",3),
                   Lon=c(145.380,145.270,135.461),
                   Lat=c(-16.4800,-5.2500,-16.0810))
  present.points<-pseudo.present.points(data,10,2,1,en.vir)

  expect_warning(niche.Model.Build(prese=present.points,absen=NULL,
                                   prese.env=NULL,absen.env=NULL,
                                   model="RF",
                                   en.vir=en.vir,bak.vir=bak.vir),
                 "The response has five or fewer unique values.  Are you sure you want to do regression?")
})

test_that("format testing on the built model2",{
  data(en.vir)
  data(bak.vir)
  data<-data.frame(species=rep("Acosmeryx anceus",3),
                   Lon=c(145.380,145.270,135.461),
                   Lat=c(-16.4800,-5.2500,-16.0810))
  present.points<-pseudo.present.points(data,10,2,1,en.vir)
  prese.env<-raster::extract(en.vir,present.points[,2:3])
  prese.env<-as.data.frame(prese.env)
  NMB.out<-niche.Model.Build(prese=NULL,absen=NULL,
                             prese.env=prese.env,absen.env=NULL,
                             model="MAXENT",
                             en.vir=en.vir,bak.vir=bak.vir)

  expect_equal(class(NMB.out$model)[1],"MaxEnt")
})

test_that("tests for abnormal conditions1",{
  data(en.vir)
  data(bak.vir)
  data<-data.frame(species=rep("Acosmeryx anceus",3),
                   Lon=c(145.380,145.270,135.461),
                   Lat=c(-16.4800,-5.2500,-16.0810))
  present.points<-data[1:2,]
  prese.env<-raster::extract(en.vir,present.points[,2:3])
  prese.env<-as.data.frame(prese.env)

  expect_warning(niche.Model.Build(prese=NULL,absen=NULL,
                                   prese.env=prese.env,absen.env=NULL,
                                   model="RF",
                                   en.vir=en.vir,bak.vir=bak.vir),
                 "prese.env has less than 3 records!\n")
})

test_that("tests for abnormal conditions2",{
  data(en.vir)
  data(bak.vir)
  data<-data.frame(species=rep("Acosmeryx anceus",3),
                   Lon=c(145.380,145.270,135.461),
                   Lat=c(-16.4800,-5.2500,-16.0810))
  present.points<-data[,1:2]

  expect_error(niche.Model.Build(prese=present.points,absen=NULL,
                                 prese.env=NULL,absen.env=NULL,
                                 model="RF",
                                 en.vir=en.vir,bak.vir=bak.vir),
               "The present data must be a dataframe with three columns")
})


test_that("tests for abnormal conditions3",{
  data(en.vir)
  data(bak.vir)
  data<-data.frame(species=rep("Acosmeryx anceus",3),
                   Lon=c(145.380,145.270,135.461),
                   Lat=c(-16.4800,-5.2500,-16.0810))
  present.points<-data
  expect_warning(niche.Model.Build(prese=present.points,absen=NULL,
                                   prese.env=NULL,absen.env=NULL,
                                   model="RF",
                                   en.vir=en.vir,bak.vir=bak.vir),
                 "The response has five or fewer unique values.  Are you sure you want to do regression?")
})

test_that("tests for abnormal conditions4",{
  data(en.vir)
  data(bak.vir)
  data<-data.frame(species=rep("Acosmeryx anceus",3),
                   Lon=c(145.380,145.270,135.461),
                   Lat=c(-16.4800,-5.2500,-16.0810))
  present.points<-data[1:2,]
  prese.env<-raster::extract(en.vir,present.points[,2:3])
  prese.env<-as.data.frame(prese.env)

  expect_warning(niche.Model.Build(prese=present.points,absen=NULL,
                                   prese.env=prese.env,absen.env=NULL,
                                   model="RF",
                                   en.vir=en.vir,bak.vir=bak.vir),
               "prese.env has less than 3 records!\n")
})

test_that("tests for abnormal conditions5",{
  data(en.vir)
  data(bak.vir)
  expect_warning(niche.Model.Build(prese=NULL,absen=NULL,
                                  prese.env=bak.vir,absen.env=NULL,
                                  model="RF",
                                  en.vir=en.vir,bak.vir=bak.vir),
                "The response has five or fewer unique values.  Are you sure you want to do regression?")
})

test_that("tests for abnormal conditions6",{
  data(en.vir)
  data(bak.vir)
  expect_warning(niche.Model.Build(prese=NULL,absen=NULL,
                                   prese.env=bak.vir,absen.env=bak.vir,
                                   model="RF",
                                   en.vir=en.vir,bak.vir=bak.vir),
                 "The response has five or fewer unique values.  Are you sure you want to do regression?")
})

test_that("tests for abnormal conditions7",{
  data(en.vir)
  data(bak.vir)
  data<-data.frame(species=rep("Acosmeryx anceus",3),
                   Lon=c(145.380,145.270,135.461),
                   Lat=c(-16.4800,-5.2500,-16.0810))
  absent.points<-pseudo.absent.points(data,en.vir=en.vir,outputNum=100)
  expect_warning(niche.Model.Build(prese=data,absen=absent.points$lonlat,
                                   prese.env=NULL,absen.env=NULL,
                                   model="RF",
                                   en.vir=en.vir,bak.vir=bak.vir),
                 "The response has five or fewer unique values.  Are you sure you want to do regression?")
})

