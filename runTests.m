%runTests
% Runs whole-cell tests and logs results in JUnit-style XML.
%
% Example:
% >> setWarnings();
% >> setPath();
% >> generateTestFixtures();
% >> runTests();
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Stanford University
% Created: 3/3/2013
function runTests(isMakeFixtures)

%import classes
import wholecell.test.XMLTestRunDisplay;
import wholecell.test.runtests;
import wholecell.util.mkdir;

%generate fixtures
if nargin < 1 || isMakeFixtures
    generateTestFixtures();
end

%create output directory
mkdir('out/test');

%instantiate XML run monitor
monitor = XMLTestRunDisplay('Whole-cell tests', 'Whole cell simulation tests', 'out/test/results.xml');

%run and log tests
runtests(monitor, {
    'wholecell.kb'
    'wholecell.sim'
    'wholecell.util'
    });
