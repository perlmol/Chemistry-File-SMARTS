use 5.006001;
use ExtUtils::MakeMaker;

WriteMakefile(
    'NAME'		=> 'Chemistry::File::SMARTS',
    'VERSION_FROM'	=> 'SMARTS.pm', 
    'PREREQ_PM'		=> {
        'Chemistry::Mol' => '0.10',
        'Chemistry::Pattern' => '0.01',
        'Parse::RecDescent' => '0',
    }, 
);
