module.exports = function(grunt) {
    // tasks
    grunt.loadTasks('build/tasks');
    grunt.loadNpmTasks('grunt-contrib-jshint');
    // configuration
    grunt.initConfig({
        pkg: '<json:package.json>',
        meta: {
            banner: '/**\n' +
                ' * <%= pkg.description %>\n' +
                ' * @version v<%= pkg.version %> - ' + '<%= grunt.template.today("yyyy-mm-dd") %>\n' +
                ' * @link <%= pkg.homepage %>\n' +
                ' * @license MIT License, http://www.opensource.org/licenses/MIT\n' +
                ' */'
        },
        jshint: {
            files: ['Converter.js']
        }
    });

    grunt.registerTask('default', ['jshint']);
};
