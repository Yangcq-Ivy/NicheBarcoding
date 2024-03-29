% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/niche.PCA.R
\name{niche.PCA}
\alias{niche.PCA}
\title{Principal component analysis of ecological niche among unknown species and
the potential species to which they may belong}
\usage{
niche.PCA(ref.lonlat, que.lonlat, en.vir = NULL)
}
\arguments{
\item{ref.lonlat}{Data frame, longitude and latitude of the known species.}

\item{que.lonlat}{Data frame, longitude and latitude of unknown species.}

\item{en.vir}{RasterBrick, the globle bioclimate data obtained from
"raster::getData" function.}
}
\value{
A list containing inportance and loadings of the components.

A dataframe of points that within the 95\% confidence interval of
the reference dataset ecological space.

A figure shows whether the query points (blue solid circles) are
located in the 95\%CI range of the niche space of reference species.
}
\description{
Determine whether unknown species belong to a known species
through principal component analysis of ecological niche according to their
distribution information.
}
\examples{
data(en.vir)
#envir<-raster::getData("worldclim",download=FALSE,var="bio",res=2.5)
#en.vir<-raster::brick(envir)

data(LappetMoths)
ref.infor<-LappetMoths$ref.infor
que.infor<-LappetMoths$que.infor

#windows() # open a new plot window when the image format is abnormal
nPCA<-niche.PCA(ref.lonlat=ref.infor[,3:5],
                que.lonlat=que.infor[,c(2,4:5)],
                en.vir=en.vir)
nPCA$summary
nPCA$que.CI


data<-data.frame(species=rep("Acosmeryx anceus",3),
                 Lon=c(145.380,145.270,135.461),
                 Lat=c(-16.4800,-5.2500,-16.0810))
simuSites<-pseudo.present.points(data,500,4,2,en.vir)
ref.lonlat<-simuSites[1:480,]
que.lonlat<-simuSites[481:500,]

#windows() # open a new plot window when the image format is abnormal
nPCA2<-niche.PCA(ref.lonlat,que.lonlat,en.vir=en.vir)
nPCA2$summary
nPCA2$que.CI
}
\author{
Cai-qing YANG (Email: yangcq_ivy(at)163.com) and Ai-bing ZHANG
(Email:zhangab2008(at)cnu.edu.cn), Capital Normal University (CNU), Beijing,
CHINA.
}
\keyword{niche.PCA}
