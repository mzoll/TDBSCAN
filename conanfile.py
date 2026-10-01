# This file is managed by Conan, contents will be overwritten.
# To keep your changes, remove these comment lines, but the plugin won't be able to modify your requirements

from conan import ConanFile
from conan.tools.cmake import cmake_layout, CMakeToolchain

class ConanApplication(ConanFile):
    package_type = "application"
    settings = "os", "compiler", "build_type", "arch"
    generators = "CMakeDeps"

    exports_sources = "CMakeLists.txt", "cmake/*", "src/*", "test/*"

    def build_requirements(self):
        self.tool_requires("cmake/[>=3.27.0]")
        self.tool_requires("ccache/[>=4.8.3]")
        self.tool_requires("cppcheck/[>=2.12.1]")
        self.test_requires("gtest/[>=1.14.0]")

    def validate(self):
        check_min_cppstd(self, "17")

    def layout(self):
        cmake_layout(self)

    def generate(self):
        tc = CMakeToolchain(self)
        tc.user_presets_path = False
        tc.generate()

    def requirements(self):
        requirements = self.conan_data.get('requirements', [])
        for requirement in requirements:
            self.requires(requirement)