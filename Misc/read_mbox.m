function metadata = read_mbox(mboxFile, opts)
%READ_MBOX Extract metadata from a Google Takeout MBOX file.
%
% metadata = read_mbox("gmail.mbox")
%
% The MBOX is divided into byte ranges and parsed in parallel. Each worker
% opens the file independently and returns messages whose starting offsets
% fall within its assigned range.
%
% Options:
%
%   UseParallel
%       Use parallel workers.
%       Default: true
%
%   NumWorkers
%       Requested number of workers.
%       Default: 4
%
%   ChunkSizeMB
%       Amount read at once by each worker.
%       Default: 64
%
%   OverlapKB
%       Extra data read around each worker's nominal range.
%       Default: 64
%
%   MaxHeaderMB
%       Maximum accepted size of one message header.
%       Default: 0.125 (128 KB)
%
%   GrowBy
%       Number of local result rows allocated at a time.
%       Default: 20000
%
%   ProgressEverySeconds
%       Progress-reporting interval.
%       Default: 5
%
%   ParseDates
%       Convert Date headers to datetime values.
%       Default: false
%
%   DecodeSubjects
%       Attempt RFC 2047 subject decoding using JavaMail.
%       Default: false
%
% Example:
%
% metadata = read_mbox( ...
%     "C:\Temp\gmail.mbox", ...
%     NumWorkers=4, ...
%     ChunkSizeMB=64, ...
%     MaxHeaderMB=0.125, ...
%     ParseDates=true);

    arguments
        mboxFile (1,1) string

        opts.UseParallel (1,1) logical = true

        opts.NumWorkers (1,1) double ...
            {mustBePositive, mustBeInteger} = 4

        opts.ChunkSizeMB (1,1) double ...
            {mustBePositive, mustBeInteger} = 64

        opts.OverlapKB (1,1) double ...
            {mustBePositive, mustBeInteger} = 64

        opts.MaxHeaderMB (1,1) double ...
            {mustBePositive} = 0.125

        opts.GrowBy (1,1) double ...
            {mustBePositive, mustBeInteger} = 20000

        opts.ProgressEverySeconds (1,1) double ...
            {mustBePositive} = 5

        opts.ParseDates (1,1) logical = false
        opts.DecodeSubjects (1,1) logical = false
    end

    filename = char(mboxFile);
    fileInfo = dir(filename);

    if isempty(fileInfo)
        error( ...
            'read_mbox:FileNotFound', ...
            'MBOX file not found: %s', ...
            filename);
    end

    if fileInfo.isdir
        error( ...
            'read_mbox:NotAFile', ...
            'The specified path is a directory: %s', ...
            filename);
    end

    fileSize = double(fileInfo.bytes);

    if fileSize == 0
        metadata = emptyMetadataTable();
        return
    end

    settings = struct();
    settings.ChunkBytes = round(opts.ChunkSizeMB * 1024^2);
    settings.OverlapBytes = round(opts.OverlapKB * 1024);
    settings.MaxHeaderBytes = round(opts.MaxHeaderMB * 1024^2);
    settings.GrowBy = opts.GrowBy;
    settings.ProgressEverySeconds = opts.ProgressEverySeconds;
    settings.ParseDates = opts.ParseDates;
    settings.DecodeSubjects = opts.DecodeSubjects;

    fprintf( ...
        'Reading %.2f GB from:\n%s\n\n', ...
        fileSize / 1024^3, ...
        filename);

    %% Start or reuse a parallel pool

    useParallel = opts.UseParallel && opts.NumWorkers > 1;
    pool = [];

    if useParallel
        try
            pool = gcp('nocreate');

            if isempty(pool)
                try
                    pool = parpool('Threads', opts.NumWorkers);
                catch threadPoolError
                    warning( ...
                        'read_mbox:ThreadPoolUnavailable', ...
                        ['Could not start a thread pool: %s\n' ...
                         'Trying a process-based local pool instead.'], ...
                        threadPoolError.message);

                    pool = parpool('local', opts.NumWorkers);
                end
            end

        catch poolError
            warning( ...
                'read_mbox:ParallelUnavailable', ...
                ['Parallel execution could not be started: %s\n' ...
                 'Continuing in serial mode.'], ...
                poolError.message);

            useParallel = false;
            pool = [];
        end
    end

    if useParallel
        numberOfWorkers = pool.NumWorkers;

        if numberOfWorkers ~= opts.NumWorkers
            warning( ...
                'read_mbox:ExistingPoolSize', ...
                ['The existing pool has %d workers rather than the ' ...
                 'requested %d. The existing pool will be used.'], ...
                numberOfWorkers, ...
                opts.NumWorkers);
        end
    else
        numberOfWorkers = 1;
    end

    %% Parse the file

    overallTimer = tic;

    if useParallel
        fprintf( ...
            'Using %d parallel workers.\n\n', ...
            numberOfWorkers);

        rangeEdges = round( ...
            linspace(0, fileSize, numberOfWorkers + 1));

        parts = cell(numberOfWorkers, 1);

        progressQueue = parallel.pool.DataQueue;
        progressBytes = zeros(numberOfWorkers, 1);
        progressMessages = zeros(numberOfWorkers, 1);
        lastProgressPrint = 0;

        afterEach(progressQueue, @updateParallelProgress);

        parfor workerIndex = 1:numberOfWorkers
            rangeStart = rangeEdges(workerIndex);
            rangeEnd = rangeEdges(workerIndex + 1);

            parts{workerIndex} = parseMboxRange( ...
                filename, ...
                rangeStart, ...
                rangeEnd, ...
                fileSize, ...
                settings, ...
                progressQueue, ...
                workerIndex);
        end

        fprintf('\n');

        metadata = vertcat(parts{:});

    else
        fprintf('Using serial parsing.\n\n');

        metadata = parseMboxRange( ...
            filename, ...
            0, ...
            fileSize, ...
            fileSize, ...
            settings, ...
            [], ...
            1);

        fprintf('\n');
    end

    %% Sort, deduplicate and renumber

    if height(metadata) > 0
        metadata = sortrows(metadata, 'FileOffset');

        [~, uniqueRows] = unique( ...
            metadata.FileOffset, ...
            'stable');

        if numel(uniqueRows) < height(metadata)
            numberRemoved = height(metadata) - numel(uniqueRows);

            warning( ...
                'read_mbox:DuplicateOffsets', ...
                'Removed %d duplicate message offsets.', ...
                numberRemoved);

            metadata = metadata(uniqueRows, :);
        end

        metadata.MessageNumber = ...
            uint32((1:height(metadata)).');
    end

    elapsed = toc(overallTimer);

    fprintf( ...
        'Parsed %d messages in %s.\n', ...
        height(metadata), ...
        formatDuration(elapsed));

    %% Parallel progress callback

    function updateParallelProgress(data)
        workerIndex = data(1);
        workerBytesProcessed = data(2);
        workerMessageCount = data(3);

        progressBytes(workerIndex) = max( ...
            progressBytes(workerIndex), ...
            workerBytesProcessed);

        progressMessages(workerIndex) = max( ...
            progressMessages(workerIndex), ...
            workerMessageCount);

        elapsedNow = toc(overallTimer);

        if elapsedNow - lastProgressPrint < ...
                opts.ProgressEverySeconds
            return
        end

        totalBytesProcessed = min( ...
            sum(progressBytes), ...
            fileSize);

        fractionComplete = totalBytesProcessed / fileSize;

        if elapsedNow > 0 && totalBytesProcessed > 0
            aggregateRate = totalBytesProcessed / elapsedNow;
            remainingSeconds = ...
                (fileSize - totalBytesProcessed) / aggregateRate;
        else
            aggregateRate = NaN;
            remainingSeconds = NaN;
        end

        fprintf( ...
            ['\r%6.2f%% | %d messages | %6.1f MB/s | ' ...
             'elapsed %s | remaining %s'], ...
            100 * fractionComplete, ...
            round(sum(progressMessages)), ...
            aggregateRate / 1024^2, ...
            formatDuration(elapsedNow), ...
            formatDuration(remainingSeconds));

        drawnow limitrate

        lastProgressPrint = elapsedNow;
    end
end


function metadata = parseMboxRange( ...
        filename, ...
        rangeStart, ...
        rangeEnd, ...
        fileSize, ...
        settings, ...
        progressQueue, ...
        workerIndex)
%PARSEMBOXRANGE Parse messages starting within one assigned byte range.

    fid = fopen(filename, 'rb');

    if fid < 0
        error( ...
            'read_mbox:WorkerCannotOpenFile', ...
            'A worker could not open: %s', ...
            filename);
    end

    cleanupObject = onCleanup(@() fclose(fid)); %#ok<NASGU>

    overlapBytes = settings.OverlapBytes;
    maxHeaderBytes = settings.MaxHeaderBytes;
    chunkBytes = settings.ChunkBytes;

    readStart = max(0, rangeStart - overlapBytes);

    readLimit = min( ...
        fileSize, ...
        rangeEnd + overlapBytes + maxHeaderBytes);

    status = fseek(fid, readStart, 'bof');

    if status ~= 0
        error( ...
            'read_mbox:WorkerSeekFailed', ...
            'A worker could not seek to byte offset %.0f.', ...
            readStart);
    end

    %% Local output storage

    capacity = settings.GrowBy;
    numberOfMessages = 0;
    skippedHeaders = 0;

    MessageNumber = zeros(capacity, 1, 'uint32');
    FromEnvelope = strings(capacity, 1);
    From = strings(capacity, 1);
    To = strings(capacity, 1);
    Cc = strings(capacity, 1);
    Bcc = strings(capacity, 1);
    Subject = strings(capacity, 1);
    DateText = strings(capacity, 1);
    DateTime = NaT(capacity, 1, 'TimeZone', 'UTC');
    MessageID = strings(capacity, 1);
    Labels = strings(capacity, 1);
    ThreadID = strings(capacity, 1);
    FileOffset = zeros(capacity, 1, 'uint64');

    %% Parser state

    overlap = zeros(1, 0, 'uint8');
    lastSeparatorOffset = -Inf;

    pendingActive = false;
    pendingOffset = 0;
    pendingEnvelope = "";
    pendingBytes = zeros(1, 0, 'uint8');

    stopScanning = false;

    workerTimer = tic;
    nextProgressTime = settings.ProgressEverySeconds;

    while double(ftell(fid)) < readLimit
        chunkStartOffset = double(ftell(fid));

        numberToRead = min( ...
            chunkBytes, ...
            readLimit - chunkStartOffset);

        chunk = fread( ...
            fid, ...
            numberToRead, ...
            '*uint8').';

        if isempty(chunk)
            break
        end

        %% Complete a header begun in the previous chunk

        if pendingActive
            remainingCapacity = ...
                maxHeaderBytes - numel(pendingBytes);

            if remainingCapacity <= 0
                skippedHeaders = skippedHeaders + 1;
                pendingActive = false;
                pendingBytes = zeros(1, 0, 'uint8');

            else
                bytesToTake = min( ...
                    remainingCapacity, ...
                    numel(chunk));

                combinedHeader = [ ...
                    pendingBytes, ...
                    chunk(1:bytesToTake)];

                [headerBytes, headerComplete] = ...
                    extractHeaderBlock(combinedHeader);

                if headerComplete
                    appendRecord( ...
                        headerBytes, ...
                        pendingEnvelope, ...
                        pendingOffset);

                    pendingActive = false;
                    pendingBytes = zeros(1, 0, 'uint8');

                elseif numel(combinedHeader) >= ...
                        maxHeaderBytes

                    skippedHeaders = skippedHeaders + 1;
                    pendingActive = false;
                    pendingBytes = zeros(1, 0, 'uint8');

                else
                    pendingBytes = combinedHeader;
                end
            end
        end

        if pendingActive
            updateOverlap();
            reportWorkerProgress(false);
            continue
        end

        %% Locate candidate MBOX separator lines

        scanBuffer = [overlap, chunk];

        scanStartOffset = ...
            chunkStartOffset - numel(overlap);

        lfPositions = find(scanBuffer == 10);

        validLF = lfPositions( ...
            lfPositions + 5 <= numel(scanBuffer));

        if isempty(validLF)
            separatorStarts = zeros(1, 0);
        else
            isFromStart = ...
                scanBuffer(validLF + 1) == uint8('F') & ...
                scanBuffer(validLF + 2) == uint8('r') & ...
                scanBuffer(validLF + 3) == uint8('o') & ...
                scanBuffer(validLF + 4) == uint8('m') & ...
                scanBuffer(validLF + 5) == uint8(' ');

            separatorStarts = validLF(isFromStart) + 1;
        end

        % The first message may begin at byte zero rather than after LF.
        if scanStartOffset == 0 && ...
                numel(scanBuffer) >= 5 && ...
                isequal(scanBuffer(1:5), uint8('From '))

            separatorStarts = [1, separatorStarts];
        end

        separatorStarts = unique( ...
            separatorStarts, ...
            'stable');

        for separatorStart = separatorStarts
            absoluteOffset = ...
                scanStartOffset + separatorStart - 1;

            if absoluteOffset <= lastSeparatorOffset
                continue
            end

            lineEndIndex = find( ...
                lfPositions >= separatorStart, ...
                1, ...
                'first');

            if isempty(lineEndIndex)
                continue
            end

            lineEnd = lfPositions(lineEndIndex);

            separatorBytes = ...
                scanBuffer(separatorStart:lineEnd - 1);

            if ~isempty(separatorBytes) && ...
                    separatorBytes(end) == 13

                separatorBytes(end) = [];
            end

            if ~isLikelyMboxSeparator(separatorBytes)
                continue
            end

            lastSeparatorOffset = absoluteOffset;

            if absoluteOffset >= rangeEnd
                stopScanning = true;
                break
            end

            if absoluteOffset < rangeStart
                continue
            end

            if numel(separatorBytes) > 5
                envelope = string( ...
                    strtrim(char(separatorBytes(6:end))));
            else
                envelope = "";
            end

            headerStart = lineEnd + 1;

            if headerStart > numel(scanBuffer)
                pendingActive = true;
                pendingOffset = absoluteOffset;
                pendingEnvelope = envelope;
                pendingBytes = zeros(1, 0, 'uint8');
                break
            end

            availableBytes = scanBuffer(headerStart:end);

            bytesToInspect = min( ...
                numel(availableBytes), ...
                maxHeaderBytes);

            [headerBytes, headerComplete] = ...
                extractHeaderBlock( ...
                    availableBytes(1:bytesToInspect));

            if headerComplete
                appendRecord( ...
                    headerBytes, ...
                    envelope, ...
                    absoluteOffset);

            elseif numel(availableBytes) >= ...
                    maxHeaderBytes

                skippedHeaders = skippedHeaders + 1;

            else
                pendingActive = true;
                pendingOffset = absoluteOffset;
                pendingEnvelope = envelope;
                pendingBytes = availableBytes;
                break
            end
        end

        updateOverlap();
        reportWorkerProgress(false);

        if stopScanning && ~pendingActive
            break
        end
    end

    if pendingActive
        skippedHeaders = skippedHeaders + 1;
    end

    reportWorkerProgress(true);

    %% Construct the local table

    if numberOfMessages == 0
        metadata = emptyMetadataTable();

    else
        rows = 1:numberOfMessages;

        metadata = table( ...
            MessageNumber(rows), ...
            FromEnvelope(rows), ...
            From(rows), ...
            To(rows), ...
            Cc(rows), ...
            Bcc(rows), ...
            Subject(rows), ...
            DateText(rows), ...
            DateTime(rows), ...
            MessageID(rows), ...
            Labels(rows), ...
            ThreadID(rows), ...
            FileOffset(rows));

        metadata.Properties.VariableNames = ...
            metadataVariableNames();
    end

    if skippedHeaders > 0
        warning( ...
            'read_mbox:SkippedHeaders', ...
            ['Worker %d skipped %d malformed or oversized ' ...
             'message headers.'], ...
            workerIndex, ...
            skippedHeaders);
    end

    %% Nested worker functions

    function appendRecord(headerBytes, envelope, offset)
        record = parseHeaderBlock(headerBytes);

        if settings.DecodeSubjects
            record.Subject = ...
                decodeMimeHeader(record.Subject);
        end

        numberOfMessages = numberOfMessages + 1;

        if numberOfMessages > capacity
            capacity = capacity + settings.GrowBy;

            MessageNumber(capacity, 1) = uint32(0);
            FromEnvelope(capacity, 1) = "";
            From(capacity, 1) = "";
            To(capacity, 1) = "";
            Cc(capacity, 1) = "";
            Bcc(capacity, 1) = "";
            Subject(capacity, 1) = "";
            DateText(capacity, 1) = "";
            DateTime(capacity, 1) = NaT('TimeZone', 'UTC');
            MessageID(capacity, 1) = "";
            Labels(capacity, 1) = "";
            ThreadID(capacity, 1) = "";
            FileOffset(capacity, 1) = uint64(0);
        end

        MessageNumber(numberOfMessages) = ...
            uint32(numberOfMessages);

        FromEnvelope(numberOfMessages) = envelope;
        From(numberOfMessages) = record.From;
        To(numberOfMessages) = record.To;
        Cc(numberOfMessages) = record.Cc;
        Bcc(numberOfMessages) = record.Bcc;
        Subject(numberOfMessages) = record.Subject;
        DateText(numberOfMessages) = record.DateText;
        MessageID(numberOfMessages) = record.MessageID;
        Labels(numberOfMessages) = record.Labels;
        ThreadID(numberOfMessages) = record.ThreadID;
        FileOffset(numberOfMessages) = uint64(offset);

        if settings.ParseDates
            DateTime(numberOfMessages) = ...
                parseEmailDate(record.DateText);
        end
    end

    function updateOverlap()
        numberToRetain = min( ...
            overlapBytes, ...
            numel(chunk));

        if numberToRetain == 0
            overlap = zeros(1, 0, 'uint8');
        else
            overlap = ...
                chunk(end - numberToRetain + 1:end);
        end
    end

    function reportWorkerProgress(forceReport)
        elapsedWorker = toc(workerTimer);

        if ~forceReport && ...
                elapsedWorker < nextProgressTime
            return
        end

        if forceReport
            nominalBytesProcessed = ...
                rangeEnd - rangeStart;
        else
            currentPosition = double(ftell(fid));

            nominalBytesProcessed = max( ...
                0, ...
                min(currentPosition, rangeEnd) - rangeStart);
        end

        if isempty(progressQueue)
            fractionComplete = ...
                nominalBytesProcessed / ...
                max(1, rangeEnd - rangeStart);

            if elapsedWorker > 0 && ...
                    nominalBytesProcessed > 0

                readRate = ...
                    nominalBytesProcessed / elapsedWorker;

                remainingSeconds = ...
                    ((rangeEnd - rangeStart) - ...
                     nominalBytesProcessed) / readRate;
            else
                readRate = NaN;
                remainingSeconds = NaN;
            end

            fprintf( ...
                ['\r%6.2f%% | %d messages | %6.1f MB/s | ' ...
                 'elapsed %s | remaining %s'], ...
                100 * fractionComplete, ...
                numberOfMessages, ...
                readRate / 1024^2, ...
                formatDuration(elapsedWorker), ...
                formatDuration(remainingSeconds));

            drawnow limitrate

        else
            send( ...
                progressQueue, ...
                [workerIndex, nominalBytesProcessed, numberOfMessages]);
        end

        nextProgressTime = ...
            elapsedWorker + settings.ProgressEverySeconds;
    end
end


function tf = isLikelyMboxSeparator(lineBytes)
%ISLIKELYMBOXSEPARATOR Test whether a line resembles an MBOX separator.

    if numel(lineBytes) < 10 || ...
            ~isequal(lineBytes(1:5), uint8('From '))

        tf = false;
        return
    end

    line = char(lineBytes);

    expression = ...
        ['^From\s+\S+\s+' ...
         '(Mon|Tue|Wed|Thu|Fri|Sat|Sun)\s+' ...
         '(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)\s+' ...
         '\d{1,2}\s+' ...
         '\d{2}:\d{2}(?::\d{2})?\s+' ...
         '(?:\S+\s+)?' ...
         '\d{4}\s*$'];

    tf = ~isempty( ...
        regexp(line, expression, 'once'));
end


function [headerBytes, complete] = extractHeaderBlock(bytes)
%EXTRACTHEADERBLOCK Return bytes before the first blank header line.
%
% Supports LF LF and CR LF CR LF.

    headerBytes = bytes;
    complete = false;

    numberOfBytes = numel(bytes);

    if numberOfBytes < 2
        return
    end

    lfLfPosition = find( ...
        bytes(1:end - 1) == 10 & ...
        bytes(2:end) == 10, ...
        1, ...
        'first');

    if numberOfBytes >= 4
        crLfCrLfPosition = find( ...
            bytes(1:end - 3) == 13 & ...
            bytes(2:end - 2) == 10 & ...
            bytes(3:end - 1) == 13 & ...
            bytes(4:end) == 10, ...
            1, ...
            'first');
    else
        crLfCrLfPosition = [];
    end

    if isempty(lfLfPosition) && ...
            isempty(crLfCrLfPosition)
        return
    end

    if isempty(lfLfPosition)
        terminatorStart = crLfCrLfPosition;
    elseif isempty(crLfCrLfPosition)
        terminatorStart = lfLfPosition;
    else
        terminatorStart = min( ...
            lfLfPosition, ...
            crLfCrLfPosition);
    end

    headerBytes = bytes(1:terminatorStart - 1);
    complete = true;
end


function record = parseHeaderBlock(headerBytes)
%PARSEHEADERBLOCK Parse selected RFC and Gmail-specific headers.

    record = struct( ...
        'From', "", ...
        'To', "", ...
        'Cc', "", ...
        'Bcc', "", ...
        'Subject', "", ...
        'DateText', "", ...
        'MessageID', "", ...
        'Labels', "", ...
        'ThreadID', "");

    if isempty(headerBytes)
        return
    end

    headerText = char(headerBytes);

    headerText = strrep( ...
        headerText, ...
        sprintf('\r\n'), ...
        sprintf('\n'));

    headerText = strrep( ...
        headerText, ...
        sprintf('\r'), ...
        sprintf('\n'));

    rawLines = regexp( ...
        headerText, ...
        '\n', ...
        'split');

    %% Unfold continuation lines

    unfoldedLines = cell(0, 1);

    for lineIndex = 1:numel(rawLines)
        line = rawLines{lineIndex};

        isContinuation = ...
            ~isempty(line) && ...
            (line(1) == ' ' || line(1) == char(9));

        if isContinuation && ~isempty(unfoldedLines)
            unfoldedLines{end} = [ ...
                unfoldedLines{end}, ...
                ' ', ...
                strtrim(line)];
        else
            unfoldedLines{end + 1, 1} = line; %#ok<AGROW>
        end
    end

    %% Extract selected fields

    for lineIndex = 1:numel(unfoldedLines)
        line = unfoldedLines{lineIndex};

        colonPosition = find( ...
            line == ':', ...
            1, ...
            'first');

        if isempty(colonPosition)
            continue
        end

        headerName = lower( ...
            strtrim(line(1:colonPosition - 1)));

        headerValue = string( ...
            strtrim(line(colonPosition + 1:end)));

        switch headerName
            case 'from'
                record.From = appendHeaderValue( ...
                    record.From, ...
                    headerValue);

            case 'to'
                record.To = appendHeaderValue( ...
                    record.To, ...
                    headerValue);

            case 'cc'
                record.Cc = appendHeaderValue( ...
                    record.Cc, ...
                    headerValue);

            case 'bcc'
                record.Bcc = appendHeaderValue( ...
                    record.Bcc, ...
                    headerValue);

            case 'subject'
                record.Subject = appendHeaderValue( ...
                    record.Subject, ...
                    headerValue);

            case 'date'
                record.DateText = appendHeaderValue( ...
                    record.DateText, ...
                    headerValue);

            case 'message-id'
                record.MessageID = appendHeaderValue( ...
                    record.MessageID, ...
                    headerValue);

            case 'x-gmail-labels'
                record.Labels = appendHeaderValue( ...
                    record.Labels, ...
                    headerValue);

            case 'x-gm-thrid'
                record.ThreadID = appendHeaderValue( ...
                    record.ThreadID, ...
                    headerValue);
        end
    end
end


function value = appendHeaderValue(existingValue, newValue)
%APPENDHEADERVALUE Combine repeated instances of the same header.

    if strlength(existingValue) == 0
        value = newValue;
    else
        value = existingValue + "; " + newValue;
    end
end


function value = decodeMimeHeader(value)
%DECODEMIMEHEADER Attempt RFC 2047 decoding using JavaMail.

    if strlength(value) == 0
        return
    end

    try
        value = string( ...
            javax.mail.internet.MimeUtility.decodeText( ...
                char(value)));
    catch
        % Retain the original value if JavaMail is unavailable.
    end
end


function dateTime = parseEmailDate(dateText)
%PARSEEMAILDATE Convert common RFC email date formats to UTC datetime.

    dateTime = NaT('TimeZone', 'UTC');

    if strlength(dateText) == 0
        return
    end

    cleaned = regexprep( ...
        char(dateText), ...
        '\s*\([^)]*\)\s*$', ...
        '');

    formats = {
        'eee, dd MMM yyyy HH:mm:ss Z'
        'eee, d MMM yyyy HH:mm:ss Z'
        'dd MMM yyyy HH:mm:ss Z'
        'd MMM yyyy HH:mm:ss Z'
        'eee, dd MMM yyyy HH:mm Z'
        'eee, d MMM yyyy HH:mm Z'
    };

    for formatIndex = 1:numel(formats)
        try
            dateTime = datetime( ...
                cleaned, ...
                'InputFormat', formats{formatIndex}, ...
                'Locale', 'en_US', ...
                'TimeZone', 'UTC');

            return
        catch
            % Try the next format.
        end
    end
end


function metadata = emptyMetadataTable()
%EMPTYMETADATATABLE Create a zero-row table with the expected types.

    metadata = table( ...
        zeros(0, 1, 'uint32'), ...
        strings(0, 1), ...
        strings(0, 1), ...
        strings(0, 1), ...
        strings(0, 1), ...
        strings(0, 1), ...
        strings(0, 1), ...
        strings(0, 1), ...
        NaT(0, 1, 'TimeZone', 'UTC'), ...
        strings(0, 1), ...
        strings(0, 1), ...
        strings(0, 1), ...
        zeros(0, 1, 'uint64'));

    metadata.Properties.VariableNames = ...
        metadataVariableNames();
end


function names = metadataVariableNames()
%METADATAVARIABLENAMES Return output names as character vectors.

    names = { ...
        'MessageNumber', ...
        'FromEnvelope', ...
        'From', ...
        'To', ...
        'Cc', ...
        'Bcc', ...
        'Subject', ...
        'DateText', ...
        'DateTime', ...
        'MessageID', ...
        'Labels', ...
        'ThreadID', ...
        'FileOffset'};
end


function text = formatDuration(seconds)
%FORMATDURATION Format seconds as HH:MM:SS.

    if ~isfinite(seconds) || seconds < 0
        text = '--:--:--';
        return
    end

    seconds = round(seconds);

    hours = floor(seconds / 3600);
    minutes = floor(mod(seconds, 3600) / 60);
    remainingSeconds = mod(seconds, 60);

    text = sprintf( ...
        '%02d:%02d:%02d', ...
        hours, ...
        minutes, ...
        remainingSeconds);
end
