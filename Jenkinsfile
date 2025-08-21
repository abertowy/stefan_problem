pipeline {
  agent {
    dockerfile true
  }
  stages {
    stage('whoami') {
      steps {
        sh 'whoami'
        sh 'pwd'
      }
    }
  }
}