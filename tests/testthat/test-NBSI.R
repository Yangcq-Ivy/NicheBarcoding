# Test Start ###
test_that("tests of different parameter combinations",{
  data(en.vir)
  data(bak.vir)
  library(ape)
  data(LappetMoths)
  ref.seq<-LappetMoths$ref.seq
  que.seq<-LappetMoths$que.seq

  expect_warning(NBSI(ref.seq,que.seq,ref.add=NULL,
                      model="MAXENT",
                      en.vir=en.vir,bak.vir=bak.vir),
                 "essentially perfect fit: summary may be unreliable")

  ref.add<-LappetMoths$ref.add
  expect_warning(NBSI(ref.seq,que.seq,ref.add=ref.add,
                      model="MAXENT",
                      en.vir=en.vir,bak.vir=bak.vir),
                 "essentially perfect fit: summary may be unreliable")
})

test_that("tests for abnormal conditions1",{
  data(en.vir)
  data(bak.vir)
  library(ape)
  data(LappetMoths)
  ref.seq<-LappetMoths$ref.seq
  que.seq<-LappetMoths$que.seq

  expect_warning(NBSI(ref.seq,que.seq=que.seq[1,],ref.add=NULL,
                      model="RF",
                      en.vir=en.vir,bak.vir=bak.vir),
                 "The response has five or fewer unique values.  Are you sure you want to do regression?")

  expect_warning(NBSI(ref.seq,que.seq=que.seq[1,],ref.add=NULL,
                       model="RF",
                      en.vir=en.vir,bak.vir=bak.vir),
                 "The response has five or fewer unique values.  Are you sure you want to do regression?")

  expect_warning(NBSI(ref.seq,que.seq=que.seq,ref.add=NULL,
                      model="RF",
                      en.vir=en.vir,bak.vir=bak.vir),
                 "The response has five or fewer unique values.  Are you sure you want to do regression?")
})

test_that("tests for abnormal conditions2",{
  data(en.vir)
  data(bak.vir)
  library(ape)
  data(LappetMoths)
  ref.seq<-LappetMoths$ref.seq
  que.seq<-LappetMoths$que.seq
  rownames(que.seq)<-gsub("[0-9\\.\\ \\-]*$","0 0",rownames(que.seq))

  expect_error(NBSI(ref.seq,que.seq,ref.add=NULL,
                    model="MAXENT",
                    en.vir=en.vir,bak.vir=bak.vir),
               "No variables can be extracted from que.infor!")

  expect_error(NBSI(ref.seq,que.seq[1,],ref.add=NULL,
                    model="MAXENT",
                    en.vir=en.vir,bak.vir=bak.vir),
               "No variables can be extracted from que.infor!")
})


