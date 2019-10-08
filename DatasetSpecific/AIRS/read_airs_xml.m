function [MetaData,GranuleID,AirsInfo]= read_airs_xml(FilePath);

%reads in the key information from an AIRS xml granule descriptor file

AirsInfo = xml2struct(FilePath);

West  = str2num(AirsInfo.S4PAGranuleMetaDataFile.SpatialDomainContainer.HorizontalSpatialDomainContainer.BoundingRectangle.WestBoundingCoordinate.Text);
East  = str2num(AirsInfo.S4PAGranuleMetaDataFile.SpatialDomainContainer.HorizontalSpatialDomainContainer.BoundingRectangle.EastBoundingCoordinate.Text );
North = str2num(AirsInfo.S4PAGranuleMetaDataFile.SpatialDomainContainer.HorizontalSpatialDomainContainer.BoundingRectangle.NorthBoundingCoordinate.Text);
South = str2num(AirsInfo.S4PAGranuleMetaDataFile.SpatialDomainContainer.HorizontalSpatialDomainContainer.BoundingRectangle.SouthBoundingCoordinate.Text);

StartDate =  AirsInfo.S4PAGranuleMetaDataFile.RangeDateTime.RangeBeginningDate.Text;
StartTime =  AirsInfo.S4PAGranuleMetaDataFile.RangeDateTime.RangeBeginningTime.Text;
EndDate   =  AirsInfo.S4PAGranuleMetaDataFile.RangeDateTime.RangeEndingDate.Text;
EndTime   =  AirsInfo.S4PAGranuleMetaDataFile.RangeDateTime.RangeEndingTime.Text;


MinTime = datenum([StartDate,' ',StartTime(1:end-2)]);
MaxTime = datenum([  EndDate,' ',EndTime(  1:end-2)]);


GranuleID = AirsInfo.S4PAGranuleMetaDataFile.DataGranule.GranuleID.Text;
MetaData = [South,West,North,East,MinTime,MaxTime];


end