BuildRequires: fftw3-devel tcl-devel make gcc

Summary: A program to simulate polymers
Name: Espresso
Version: 2.0.2
Release: d
Source: %{name}-%{version}.tar.gz
License: Espresso Licence
Group: Applications/Engineering
BuildRoot: %{_builddir}/%{name}-root

%description
A program to simulate polymers and more. Details....

%prep
%setup -q
%build
export ESPRESSO_SOURCE=%{_builddir}/%{name}-%{version}
export ESPRESSO_SCRIPTS=$ESPRESSO_SOURCE/scripts
./configure --with-mpi=fake
make
%install
export ESPRESSO_SOURCE=%{_builddir}/%{name}-%{version}
export ESPRESSO_SCRIPTS=$ESPRESSO_SOURCE/scripts
export conf=`$ESPRESSO_SOURCE/config/config.guess`
install -D obj-$conf/Espresso_bin ${RPM_BUILD_ROOT}/usr/local/Espresso/Espresso_bin
install -D start_Espresso ${RPM_BUILD_ROOT}/usr/local/Espresso/start_Espresso
mkdir -p ${RPM_BUILD_ROOT}/usr/local/Espresso/scripts
install -D scripts/* ${RPM_BUILD_ROOT}/usr/local/Espresso/scripts
mkdir -p ${RPM_BUILD_ROOT}/usr/local/bin
ln -s ${RPM_BUILD_ROOT}/usr/local/Espresso/start_Espresso ${RPM_BUILD_ROOT}/usr/local/bin/Espresso
%clean
rm -rf ${RPM_BUILD_ROOT}
%files
%defattr(-,root,root)
/usr/local/Espresso/start_Espresso
/usr/local/Espresso/Espresso_bin
/usr/local/Espresso/scripts
%changelog
* Tue Feb 27 2007 - Christoph Junghans <junghans@mpip-mainz.mpg.de>
- initial package: version 2.0.2