import { HttpClient, HttpHeaders, HttpParams } from '@angular/common/http';
import { Injectable } from '@angular/core';

export interface Res {
    mature: string[];
    pre: string[];
}

@Injectable({
    providedIn: 'root'
})

export class RequestService {

    constructor(private http: HttpClient) { }

    trimmingStep(body: string) {
        const path = 'http://localhost:8080/trimming';
        return this.http.post(
            path,
            body,
            {
                headers: new HttpHeaders({
                    'Content-Type': 'application/json',
                }),
                responseType: 'text'
            });
    }

    shortStack(body: string) {
        const path = 'http://localhost:8080/shortstack';
        return this.http.post(
            path,
            body,
            {
                headers: new HttpHeaders({
                    'Content-Type': 'application/json',
                }),
                responseType: 'text'
            });
    }

    mirnaAnalysis(body: string) {
        const path = 'http://localhost:8080/mirna';
        return this.http.post(
            path,
            body,
            {
                headers: new HttpHeaders({
                    'Content-Type': 'application/json',
                }),
                responseType: 'text'
            });
    }

    mirnaDetection() {
        const path = 'http://localhost:8080/mirnas';
        return this.http.get<Res>(path);
    }

    mirnaInformation(path: string) {
        return this.http.get(
            path,
            {
                responseType: 'text'
            });
    }
}
