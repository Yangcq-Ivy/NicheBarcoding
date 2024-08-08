# Test Start ###
test_that("test when ref and que are NULL",{
  data(LappetMoths)
  barcode.identi.result<-LappetMoths$barcode.identi.result

  expect_error(NBSI2(ref.infor=NULL,que.infor=NULL,
                     ref.env=NULL,que.env=NULL,
                     barcode.identi.result=barcode.identi.result,
                     model="MAXENT",
                     en.vir=NULL,bak.vir=NULL),
               "Please check the input data!")
})

test_that("tests of different parameter combinations",{
  data(en.vir)
  data(bak.vir)
  data(LappetMoths)
  barcode.identi.result<-LappetMoths$barcode.identi.result

  ref.infor<-LappetMoths$ref.infor
  que.infor<-LappetMoths$que.infor
  expect_warning(NBSI2(ref.infor=ref.infor,que.infor=que.infor,
                       barcode.identi.result=barcode.identi.result,
                       model="MAXENT",
                       en.vir=en.vir,bak.vir=bak.vir),
                 "essentially perfect fit: summary may be unreliable")

  ref.env<-LappetMoths$ref.env
  que.env<-LappetMoths$que.env
  expect_warning(NBSI2(ref.env=ref.env,que.env=que.env,
                       barcode.identi.result=barcode.identi.result,
                       model="MAXENT",
                       en.vir=en.vir,bak.vir=bak.vir),
                 "essentially perfect fit: summary may be unreliable")
})

test_that("tests for abnormal conditions1",{
  data(en.vir)
  data(bak.vir)
  data(LappetMoths)
  barcode.identi.result<-LappetMoths$barcode.identi.result
  ref.infor<-LappetMoths$ref.infor
  que.infor<-LappetMoths$que.infor
  expect_warning(NBSI2(ref.infor=ref.infor,que.infor=que.infor[1,],
                       barcode.identi.result=barcode.identi.result,
                       model="MAXENT",
                       en.vir=en.vir,bak.vir=bak.vir),
                 "essentially perfect fit: summary may be unreliable")

  ref.env<-LappetMoths$ref.env
  que.env<-LappetMoths$que.env
  expect_warning(NBSI2(ref.env=ref.env,que.env=que.env[1,],
                       barcode.identi.result=barcode.identi.result,
                       model="MAXENT",
                       en.vir=en.vir,bak.vir=bak.vir),
                 "essentially perfect fit: summary may be unreliable")
})

test_that("tests for abnormal conditions2",{
  data(en.vir)
  data(bak.vir)
  data(LappetMoths)
  barcode.identi.result<-LappetMoths$barcode.identi.result

  ref.infor<-LappetMoths$ref.infor
  que.infor<-LappetMoths$que.infor
  expect_warning(NBSI2(ref.infor=ref.infor,que.infor=que.infor,
                       barcode.identi.result=barcode.identi.result,
                       model="MAXENT",
                       en.vir=en.vir,bak.vir=bak.vir),
                 "essentially perfect fit: summary may be unreliable")
})

test_that("tests for abnormal conditions3",{
  data(en.vir)
  data(bak.vir)
  data(LappetMoths)
  barcode.identi.result<-LappetMoths$barcode.identi.result
  ref.infor<-LappetMoths$ref.infor
  que.env<-LappetMoths$que.env

  expect_error(NBSI2(ref.infor=ref.infor,que.env=que.env,
                       barcode.identi.result=barcode.identi.result,
                       model="RF",
                       en.vir=en.vir,bak.vir=bak.vir),
                 "There is no matching ecological variable information between REF and QUE!")
})


