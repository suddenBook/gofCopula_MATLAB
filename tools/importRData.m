function importRData(sourceFolder, outputFolder)
%IMPORTRDATA Create MAT datasets from exact binary arrays exported by R.
%   IMPORTRDATA(SOURCEFOLDER,OUTPUTFOLDER) reads output from exportRData.R.

arguments
    sourceFolder (1,1) string {mustBeFolder}
    outputFolder (1,1) string
end
if ~isfolder(outputFolder)
    mkdir(outputFolder);
end

years = readNumericLines(fullfile(sourceFolder,"Banks_years.txt"));
rowCounts = readNumericLines(fullfile(sourceFolder,"Banks_rows.txt"));
names = readTextLines(fullfile(sourceFolder,"Banks_variables.txt")).';
data = readAnnualData(sourceFolder,"Banks",years,rowCounts,numel(names));
Banks = struct("Years",years,"VariableNames",names,"Data",{data});
Provenance = provenance("Banks.RData","Yahoo Finance");
save(fullfile(outputFolder,"Banks.mat"),"Banks","Provenance","-v7");

years = readNumericLines(fullfile(sourceFolder,"CryptoCurrencies_years.txt"));
rowCounts = readNumericLines(fullfile(sourceFolder,"CryptoCurrencies_rows.txt"));
names = readTextLines(fullfile(sourceFolder,"CryptoCurrencies_variables.txt")).';
data = readAnnualData(sourceFolder,"CryptoCurrencies",years,rowCounts,numel(names));
CryptoCurrencies = struct( ... %#ok<NASGU>
    "Years",years,"VariableNames",names,"Data",{data});
Provenance = provenance("CryptoCurrencies.RData", ... %#ok<NASGU>
    "CoinMetrics (https://coinmetrics.io/)");
save(fullfile(outputFolder,"CryptoCurrencies.mat"), ...
    "CryptoCurrencies","Provenance","-v7");

dimensions = readNumericLines(fullfile(sourceFolder,"IndexReturns2D_dim.txt"));
IndexReturns2D = readMatrix( ... %#ok<NASGU>
    fullfile(sourceFolder,"IndexReturns2D.bin"),dimensions(1),dimensions(2));
VariableNames = readTextLines( ... %#ok<NASGU>
    fullfile(sourceFolder,"IndexReturns2D_variables.txt")).';
Provenance = provenance("IndexReturns2D.RData", ... %#ok<NASGU>
    "Erste Bank AG, Vienna, Austria");
save(fullfile(outputFolder,"IndexReturns2D.mat"), ...
    "IndexReturns2D","VariableNames","Provenance","-v7");

dimensions = readNumericLines(fullfile(sourceFolder,"IndexReturns3D_dim.txt"));
IndexReturns3D = readMatrix( ... %#ok<NASGU>
    fullfile(sourceFolder,"IndexReturns3D.bin"),dimensions(1),dimensions(2));
VariableNames = readTextLines( ... %#ok<NASGU>
    fullfile(sourceFolder,"IndexReturns3D_variables.txt")).';
Provenance = provenance("IndexReturns3D.RData", ... %#ok<NASGU>
    "Erste Bank AG, Vienna, Austria");
save(fullfile(outputFolder,"IndexReturns3D.mat"), ...
    "IndexReturns3D","VariableNames","Provenance","-v7");
end

function data = readAnnualData(folder,prefix,years,rowCounts,columnCount)
data = cell(numel(years),1);
for k = 1:numel(years)
    file = fullfile(folder,prefix + "_" + string(years(k)) + ".bin");
    data{k} = readMatrix(file,rowCounts(k),columnCount);
end
end

function matrix = readMatrix(file,rowCount,columnCount)
fileIdentifier = fopen(file,"r","ieee-le");
if fileIdentifier < 0
    error("gofcopula:data:CannotOpen","Cannot open %s.",file);
end
cleanup = onCleanup(@() fclose(fileIdentifier));
values = fread(fileIdentifier,rowCount*columnCount,"double=>double");
if numel(values) ~= rowCount*columnCount
    error("gofcopula:data:UnexpectedLength", ...
        "%s does not contain the expected number of doubles.",file);
end
matrix = reshape(values,rowCount,columnCount);
clear cleanup
end

function values = readNumericLines(file)
values = str2double(readTextLines(file));
if any(isnan(values))
    error("gofcopula:data:InvalidMetadata", ...
        "%s contains nonnumeric metadata.",file);
end
end

function values = readTextLines(file)
values = readlines(file,EmptyLineRule="skip");
end

function value = provenance(originalFile,source)
value = struct( ...
    "OriginalFile",originalFile, ...
    "OriginalPackage","gofCopula 0.4-3", ...
    "Source",source, ...
    "Conversion", ...
        "R writeBin(size=8, endian=little) to MATLAB double; exact binary64 transfer");
end
