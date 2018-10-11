# Christopher Iliffe Sprague
# christopher.iliffe.sprague@gmail.com

import yaml, re, bunch, pandas

class analyser(object):

    # https://intra.kth.se/polopoly_fs/1.830206!/Datalogi%20ASP%20fastst%C3%A4lld%202016-12-13.pdf

    def __init__(self, courses):

        self.courses = [bunch.Bunch(course) for course in courses]

        self.csc_courses = [
            course
            for course
            in self.courses
            if 'DD'
            in course.name
        ]
        self.third_cycle_courses = [
            course
            for course
            in self.courses
            if int(re.split(r'(^[^\d]+)', course.code)[2][0])
            == 3
        ]
        self.first_cycle_courses = [
            course
            for course
            in self.courses
            if int(re.split(r'(^[^\d]+)', course.code)[2][0])
            == 1
        ]
        self.third_cycle_csc_courses = [
            course
            for course
            in self.third_cycle_courses
            if 'DD'
            in course.code
        ]

        self.credits = sum([
            course.credits for course in self.courses
        ])
        self.third_cycle_credits = sum([
            course.credits for course in self.third_cycle_courses
        ])
        self.first_cycle_credits = sum([
            course.credits for course in self.first_cycle_courses
        ])
        self.third_cycle_csc_credits = sum([
            course.credits for course in self.third_cycle_csc_courses
        ])




    def __str__(self):
        tables = self()
        return str(tables[0]) + "\n" + str(tables[1])

    def __call__(self):

        ind0 = ["≥60 credits", "≥45 3rd cycle credits", "≤10 1st cycle credits", "≥30 3rd cycle CSC credits"]
        data0 = [self.credits, self.third_cycle_credits, self.first_cycle_credits, self.third_cycle_csc_credits]
        data0 = pandas.DataFrame(data0, index=ind0, columns=["credits"])

        ind1 = [course.name for course in self.courses]
        col1 = [key for key in self.courses[0].keys() if key != 'name']
        data1 = pandas.DataFrame([[course[key] for key in col1] for course in self.courses], index=ind1, columns=col1)

        return (data0, data1)





if __name__ == "__main__":

    with open("courses.yaml", "r") as f:
        courses = yaml.load(f)

    report = analyser(courses)
    print(report)
