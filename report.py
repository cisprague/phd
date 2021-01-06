# Christopher Iliffe Sprague
# sprague@kth.se

import yaml

def report(fname, ladok=False, completed=False):

    # course data
    courses = yaml.load(open(fname, 'r'))

    # credit critera
    total_credits = 0
    first_cycle_credits = 0
    third_cycle_credits = 0
    third_cycle_csc_credits = 0

    # completed, but, un-ladoked
    topics = list()

    # compute course credits
    for course in courses:

        # completed, but, un-ladoked
        if not course['ladok'] and course['completed']:
            topics.append(course)

        # in Ladok and completed
        if course['ladok'] if ladok else True and course['completed'] if completed else True:

            # all courses
            total_credits += course['credits']

            # first-cycle courses
            if course['code'][2] == '1':
                first_cycle_credits += course['credits']

            # third-cycle course
            if course['code'][2] == '3':
                third_cycle_credits += course['credits']

                # third-cycle CSE courses
                if course['code'][:2] == 'DD':
                    third_cycle_csc_credits += course['credits']

    # results (https://intra.kth.se/polopoly_fs/1.830206!/Datalogi%20ASP%20fastst%C3%A4lld%202016-12-13.pdf)
    print('\nCourses')
    print('>= 60 credits: {} ({})'.format(True if total_credits >= 60 else False, total_credits))
    print('<= 10 1st cycle credits: {} ({})'.format(True if first_cycle_credits <= 10 else False, first_cycle_credits))
    print('>= 45 3rd cycle credits: {} ({})'.format(True if third_cycle_credits >= 45 else False, third_cycle_credits))
    print('>= 30 3rd cycle CSC credits: {} ({})'.format(True if third_cycle_csc_credits >= 30 else False, third_cycle_csc_credits))

    print('\nCompleted, but, un-Ladoked')
    for course in topics:
        print('{} {} ({})'.format(course['code'], course['name'], course['credits']))




if __name__ == "__main__":

    report('courses.yaml')