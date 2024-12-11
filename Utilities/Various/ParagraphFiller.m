classdef ParagraphFiller

    properties

        margin                = 4
        parlength             = 40
        preventSplitCharacter = '~' 
        
    end

    methods

        function parfill = ParagraphFiller(varargin)

            parfill = merge_options(parfill, varargin{:});
            
        end

        function lines = getLines(parfill, str)

            words = parfill.splitInWords(str);

            done = false;
            lines = {};
            while ~done
                [done, words, lines] = parfill.addOneLine(words, lines);
            end
            
        end

        function [done, words, lines] = addOneLine(parfill, words, lines)

            margin    = parfill.margin;
            parlength = parfill.parlength;
            
            if isempty(words)
                done = true;
                return
            end

            lwords = cellfun(@(word) strlength(word), words);

            lwords(1 : end - 1) = lwords(1 : end - 1) + 1; % add for one space

            clwords = cumsum(lwords);

            ind = find(clwords > parlength, 1);

            if isempty(ind)

                done = true;
                line = strjoin(words(1 : end), ' ');
                lines = {lines{:}, line};
                return
                
            elseif ind == 1

                % we have to accept this word
                line = words{ind};
                lines = {lines{:}, line};
                if numel(words) == 1
                    done = true;
                else
                    done = false;
                    words = words(ind + 1 : end);
                end
                
            else

                ind = ind - 1;

                if (clwords(ind) < parlength) && (clwords(ind + 1) <= (parlength + margin))
                    ind = ind + 1;
                end

                line = strjoin(words(1 : ind), ' ');
                lines = {lines{:}, line};

                if ind == numel(words)
                    done = true;
                else
                    done = false;
                    words = words(ind + 1 : end);
                end
                
            end
            
        end

        function words = splitInWords(parfill, str)

            words = {};

            done = false;
            while ~done
                [done, str, words] = parfill.parseOneWord(str, words);
            end
            
        end


        function [done, str, words] = parseOneWord(parfill, str, words)

            if isempty(str)
                done = true
                return
            end

            rexp = ['([ ]+|[' parfill.preventSplitCharacter '])'];
            [startIndex, endIndex] = regexp(str, rexp, 'once');
            
            if isempty(startIndex)
                words = {words{:}, str};
                done = true;
                return
            else
                done = false;
                match = str(startIndex : endIndex);
                if strcmp(match, '~')
                    str = str(endIndex + 1 : end);
                    rexp = ['[' parfill.preventSplitCharacter ']'];
                    [startIndex, endIndex] = regexp(str, rexp, 'once');
                    if startIndex > 0
                        words = {words{:}, str(1 : startIndex - 1)};
                    end
                    str  = str(endIndex + 1 : end);
                else
                    if startIndex > 1
                        words = {words{:}, str(1 : startIndex - 1)};
                    end            
                    str  = str(endIndex + 1 : end);
                end
            end
            
        end

    end

end
